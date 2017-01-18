#!/usr/bin/env python

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Sample conversion script for 1,000,000 VA data.

The 1,000,000 VA data is comprised of VCFs with every position
called.  This script converts those files to gVCF format.

Assumptions:
- one sample per input VCF
"""

import os
import re
import sys
import pdb
import argparse
import gzip

# Constants
INPUT_FILE_KEY = "map_input_file"
SAMPLE_ID_PATTERN = r"/(LP\d{7}-DNA_\w\d{2})/"

# VCF Fields
# http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
CHROM = 0
POS = 1
ID = 2
REF = 3
ALT = 4
QUAL = 5
FILTER = 6
INFO = 7
FORMAT = 8
GENOTYPE = 9


class gVCFMapper(object):
  def __init__(self, input, output, min_gq=0, min_dp=0, debug=False):
    self.input = input
    self.debug = debug
    self.min_gq = int(min_gq)
    self.min_dp = int(min_dp)
    self.output_fh = open(output, 'w')

    self.g_start_block = None
    self.g_end_block = None
    self.ref_block = None

    self.run()
    self.output_fh.close()


  def run(self):
    # Loop over each line of the VCF
    f = self.open(self.input)
    for line in f:
      line = line.rstrip()
      if not line: continue
      if line.startswith('#'):
        self.emit(line)
        continue
      fields = line.split("\t")
      if self.is_variant(fields):
        self.emit_block()
        self.emit(line)
      else:
        # Gather information about this VCF line in our non-variant region
        self.accumulate_block(fields)

    # Emit the final block, if applicable
    self.emit_block()

  def open(self, file):
    if file.endswith('.gz'):
      return gzip.open(file,'r')
    else:
      return open(file)

  def accumulate_block(self, fields, no_call=False):
    """Accumulates data from each record in a non-variant region.

    Every consecutive VCF record in a non-variant region is passed to
    this function.

    You might modify it to observe the minimum or average quality score
    for records in this region.

    You might also decide within this function that these fields should
    be part of a new VCF block record instead of the current one.  If
    so, call emit_block() and then begin accumulating state for the new
    VCF block record.

    Args:
      fields: (string array) fields from one line of the VCF

    Returns: n/a

    Side Effects:
      global variables are modified
    """

    # Check to see if the current block matches what we want to add to it
    if no_call is True and self.ref_block is True:
      self.emit_block()
    elif no_call is False and self.ref_block is False:
      self.emit_block()
    elif self.meets_filter_criteria(fields) is False:
      self.emit_block()

    # If this is the first line in a no call block set the genotype to ./.
    if no_call is True and self.ref_block is None:
      fields[GENOTYPE] = './.'
      self.ref_block = False
    elif no_call is False and self.ref_block is None:
      self.ref_block = True

    # Set start and end points of block as necessary
    if self.g_start_block is None:
      self.g_start_block = fields
      self.g_end_block = fields
    # Check if calls are adjacent
    elif not self.g_start_block[CHROM] == fields[CHROM] or \
                            self.g_end_block is not None and int(fields[POS]) > int(self.block_end_value(fields)) + 1:
      self.emit_block()
      self.g_start_block = fields
      # Emit resets ref_block, need to set again
      if no_call is True:
        self.ref_block = False
      else:
        self.ref_block = True
    else:
      self.g_end_block = fields

  def block_end_value(self, fields):
    info_dict = self.info_to_dict(fields)
    if "END" in info_dict:
      return info_dict["END"]
    else:
      return fields[POS]

  def check_values(self, existing, new, metric, min=True):
    if existing is None:
      existing = new
    else:
      type = self.metric_type(metric)
      existing = type(existing)
      new = type(new)
      if min is True:
        if existing > new:
          existing = new
      else:
        if existing < new:
          existing = new
    return existing

  def metric_type(self, metric):
    metric_key = {
      "QUAL": float,
      "DP": int,
      "MQ": float,
      "MQ0": int
    }
    if not metric in metric_key:
      raise Exception("No type exists for %s" % metric)
    return metric_key[metric]

  def emit_block(self):
    """Emits the current non-variant block, if applicable.

    The current implementation naively uses all values from the first
    VCF record in the block overwriting the INFO value with the end of
    the block.

    You might modify this function to utilize more accumulated state about
    the block such as minimum or average quality score.

    You might also merge the END INFO key/value pairs with any other
    existing INFO key/values.

    Args:
      key: (string) the key upon which to sort prior to the reduce step

    Returns
    : n/a

    Side Effects:
      global variables are modified
      a VCF line is written to stdout
    """

    if self.g_start_block is None:
      return

    block_fields = self.g_start_block

    # No need to define the end if there is only one base
    #if self.g_end_block is None:
    #  end = int(self.g_start_block[POS]) + len(self.g_start_block[REF]) - 1
    #  block_fields[INFO] = "END=" + str(end)
    #else:

    if self.g_end_block:
      end = int(self.block_end_value(self.g_end_block)) + len(self.g_end_block[REF]) - 1
      block_fields[INFO] = "END=" + str(end) #

    if self.ref_block is False:
      block_fields[FORMAT] = "GT"
      block_fields[GENOTYPE] = "./."

    line = "\t".join(block_fields)
    self.emit(line)

    # Reset our block state
    self.g_start_block = None
    self.g_end_block = None
    self.ref_block = None


  def emit(self, line):
    """Emits a VCF line to stdout.

    Args:
      key: (string) the key upon which to sort prior to the reduce step
      line: (string)

    Returns: n/a

    Side Effects:
      a VCF line is written to stdout
    """
    self.output_fh.write(line + "\n")

  def is_variant(self, fields):
    """Determines whether or not the VCF fields constitute variant.

    Args:
      fields: (string array) fields from one line of the VCF

    Returns:
      (boolean) whether or not the VCF field constitute a variant

    Side Effects:
      a VCF line is written to stdout
    """
    if "0|0" in fields[GENOTYPE] or "0/0" in fields[GENOTYPE]:
      return False
    return True

  def is_snp(self, fields):
    if len(fields[ALT]) == 1:
      return True
    elif len(fields[ALT]) > 1:
      return False
    elif len(fields[REF]) > 1:
      return False
    return None

  def is_indel(self, fields):
    if fields[ALT] == '<NON_REF>':
      return False
    if len(fields[ALT]) == 1 and len(fields[REF]) == 1:
      return False
    if len(fields[ALT]) > 1 or len(fields[REF]) > 1:
      return True

  def is_ref(self, fields):
    if fields[ALT] in ('.', '<NON_REF>'):
      return True


  def meets_filter_criteria(self, fields):
    """
     filter criteria 1: (returns true if)
     if ALT field is . (Mean a reference call) then
     check that the values if the INFO field meet the following requirements:
     1) MQ0 less than 4
     2) MQ greater than or = to 30
     3) QUAL greater than or = to 30
    """
    #if fields[ALT] != ".":
    #  if fields[FILTER] == "PASS":
    #    return True
    #  else:
    #    return False
    #### No more checking quality or reference calls ####
    #if fields[ALT] == ".":    # Means it is a reference call Now to check the values meet our threshold
    #
    #  variant_info_dict = info_to_dict(fields)
    #
    #  ## Processing the above dict to meet criteria.
    #  if float(variant_info_dict['MQ0']) < 4 and float(variant_info_dict['MQ']) >= 30 and float(fields[QUAL]) >= 30:
    #    return True
    #  else:
    #    return False
    #else:
    #  return True

    if self.is_ref(fields):
      if fields[REF] == 'N':
        return False
      call_info = self.call_info(fields)
      if "GQ" in call_info:
        if self.min_gq > int(call_info["GQ"]):
          return False
      if "DP" in call_info:
        if self.min_dp > int(call_info["DP"]):
          return False

    return True



  def info_to_dict(self, fields):
    ##Preprocessing: Takes the string value of INFO column (field[INFO]) and converts into dict as Key, Value based on the =
    info_string_items = [s for s in fields[INFO].split(';') if s]
    variant_info_dict = {}
    for item in info_string_items:
      try:
        key,value = item.split('=')
        variant_info_dict[key] = value
      except:
        pass  #Exception handling for keys without values, Specifically 'DB' in the INFO string
    return variant_info_dict

  def call_info(self, fields):
    call_dict = {}
    format = fields[FORMAT].split(":")
    values = fields[GENOTYPE].split(":")
    index = 0
    for f in format:
      call_dict[f] = values[index]
      index += 1
    return call_dict

def parse_command_line():
    parser = argparse.ArgumentParser(
        description = 'Merge reference matching blocks in a gVCF and output a new gVCF.')
    parser.add_argument("-g", "--gVCF", help="gVCF File")
    parser.add_argument("-o", "--output", help="Output File")
    parser.add_argument("--min_gq", default=0, help="Minimum GQ value for reference calls")
    parser.add_argument("--min_dp", default=0, help="Minimum DP value for reference calls")
    parser.add_argument("-d", "--debug", action='store_true', default=False,
                                help="Output debugging messages.  May be very verbose.")
    options = parser.parse_args()

    if not options.gVCF or not options.output:
      print "Missing arguments"
      parser.print_help()
      exit()

    return options

# Main
if __name__ == "__main__":
  options = parse_command_line()
  gVCFMapper(options.gVCF, options.output, options.min_gq, options.min_dp, options.debug)
