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
called.  This script filters calls and outputs a line for every position.
Filtered calls are marked with ./. rather than the original genotype.

Assumptions:
- one sample per input VCF

This is a modified version of the simple gvcf conversion script.  It uses the reduce step and a custom output format to consolidate data for each sample to a smaller number of files.  For more detail, see
    http://stackoverflow.com/questions/7449756/get-input-file-name-in-streaming-hadoop-program
    http://stackoverflow.com/questions/18541503/multiple-output-files-for-hadoop-streaming-with-python-mapper
"""

import os
import re
import sys
import pdb

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

# Using global variables in an effort to keep this script simple
g_start_block = None
g_end_block = None
ref_block = None

def main():
  """Entry point to the script."""

  global sample_id
  sample_id = None
  sample_id_re = re.compile(SAMPLE_ID_PATTERN)

  # Basic parsing of command line arguments to allow a filename
  # to be passed when running this code in the debugger.
  path = None
  file_handle = sys.stdin
  if 2 <= len(sys.argv):
    path = sys.argv[1]
    file_handle = open(path, "r")
  else:
    path = os.environ[INPUT_FILE_KEY]
    print >> sys.stderr, path

  if path is not None:
    match = sample_id_re.search(path)
    if match:
        sample_id = match.group(1)

  # Loop over each line of the VCF
  line = file_handle.readline()
  while line:
    line = line.strip()
    if not line:
      continue

    if "#" == line[0]:
      # This is a header line, skip it
      line = file_handle.readline()
      continue

    fields = line.split("\t")
    if is_variant(fields):
        if meets_filter_criteria(fields) == True:
          # This is a variant, emit the preceeding non-variant region VCF block, if
          # applicable, followed by this VCF line
          emit(sample_id, line)
        else:
          emit_no_call(fields)
    else:
      # Gather information about this VCF line in our non-variant region
      if meets_filter_criteria(fields) == True:
        emit(sample_id, line)
      else:
        emit_no_call(fields)

    line = file_handle.readline()

def emit_no_call(fields):
  fields[FORMAT] = "GT"
  fields[GENOTYPE] = "./."
  line = "\t".join(fields)
  emit(sample_id, line)

def emit(key, line):
  """Emits a VCF line to stdout.

  Args:
    key: (string) the key upon which to sort prior to the reduce step
    line: (string)

  Returns: n/a

  Side Effects:
    a VCF line is written to stdout
  """
  print "%s\t%s" % (key, line)

def is_variant(fields):
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

def meets_filter_criteria(fields):
  """
   filter criteria 1: (returns true if)
   if ALT field is . (Mean a reference call) then
   check that the values if the INFO field meet the following requirements:
   1) MQ0 less than 4
   2) MQ greater than or = to 30
   3) QUAL greater than or = to 30
  """
  if fields[ALT] != ".":
    if fields[FILTER] == "PASS":
      return True
    else:
      return False
  if fields[ALT] == ".":    # Means it is a reference call Now to check the values meet our threshold

    variant_info_dict = info_to_dict(fields)

    ## Processing the above dict to meet criteria.
    if float(variant_info_dict['MQ0']) < 4 and float(variant_info_dict['MQ']) >= 30 and float(fields[QUAL]) >= 30:
      #vcf_count =+ 1
      return True
    else:
      return False
  else:
    return True

def info_to_dict(fields):
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


if __name__ == "__main__":
  main()
