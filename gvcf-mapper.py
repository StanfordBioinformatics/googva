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

TODO(Cuiping and Denis):
 (1) Add logic for VCF block records such as depth and genotype quality.  For
     more detail, see https://sites.google.com/site/gvcftools/home/about-gvcf
     and http://www.broadinstitute.org/gatk/guide/article?id=4017
 (2) Use your favority python unit test framework to add some tests for the
     logic to collapse contiguous non-variant calls into a VCF block record.

This script can be run standalone:
   zcat chr21.vcf.gz | ./gvcf-mapper.py

Or via the debugger:
   python -mpdb ./gvcf-mapper.py chr1.vcf

It can also be run as a Hadoop Streaming job:
  hadoop jar /path/to/your/hadoop-streaming-*.jar \
  -libjars /home/deflaux/custom.jar \
  -outputformat com.custom.CustomMultiOutputFormat \
  -mapper /home/deflaux/gvcf-mapper.py -file /home/deflaux/gvcf-mapper.py \
  -reducer org.apache.hadoop.mapred.lib.IdentityReducer \
  -input inputpath \
  -output outputpath

See also https://developers.google.com/hadoop/

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

vcf_count = 0
filtered_count = 0
total_count = 0

# Using global variables in an effort to keep this script simple
g_start_block = None
g_end_block = None

def main():
  """Entry point to the script."""

  global sample_id
  sample_id = None
  sample_id_re = re.compile(SAMPLE_ID_PATTERN)

  # Basic parsing of command line arguments to allow a filename
  # to be passed when running this code in the debugger.
  file_handle = sys.stdin
  if 2 <= len(sys.argv):
    path = sys.argv[1]
    file_handle = open(path, "r")
  else:
    path = os.environ[INPUT_FILE_KEY]
    print >> sys.stderr, path

  match = sample_id_re.search(path)
  sample_id = match.group(1)

  # Loop over each line of the VCF
  line = file_handle.readline()
  while line:
    total_count =+ 1
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
          emit_block(sample_id)
          emit(sample_id, line)
    else:
      # Gather information about this VCF line in our non-variant region
      if meets_filter_criteria(fields) == True:
        accumulate_block(fields)

    line = file_handle.readline()

  # Emit the final block, if applicable
  emit_block(sample_id)

def accumulate_block(fields):
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
  global g_start_block
  global g_end_block

  if g_start_block is None:
    g_start_block = fields
  elif g_end_block is not None and int(fields[POS]) > int(g_end_block[POS]) + 1:
      emit_block(sample_id)
      g_start_block = fields
  else:
    g_end_block = fields


def emit_block(key):
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
  global g_start_block
  global g_end_block

  if g_start_block is None:
    return

  block_fields = g_start_block

  if g_end_block is None:
    block_fields[INFO] = "END=" + str(int(g_start_block[POS]) + 1)
  else:
    block_fields[INFO] = "END=" + g_end_block[POS]

  line = "\t".join(block_fields)
  emit(key, line)
  # Reset our block state
  g_start_block = None
  g_end_block = None


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
   1) MQ0 less than or = to 4
   2) MQ greater than or = to 30
   3) QUAL greater than or = to 30
  """
  if fields[ALT] != ".":
    if fields[FILTER] == "PASS":
      return True
    else:
      return False
  if fields[ALT] == ".":    # Means it is a reference call Now to check the values meet our threshold
    ##Preprocessing: Takes the string value of INFO column (field[INFO]) and converts into dict as Key, Value based on the =
    info_string_items = [s for s in fields[INFO].split(';') if s]
    variant_info_dict = {}
    for item in info_string_items:
      try:
        key,value = item.split('=')
        variant_info_dict[key] = value
      except:
        pass  #Exception handling for keys without values, Specifically 'DB' in the INFO string
    
    ## Processing the above dict to meet criteria.
    if float(variant_info_dict['MQ0']) <= 4 and float(variant_info_dict['MQ']) >= 30 and float(fields[QUAL]) >= 30:
      vcf_count =+ 1
      return True
    else:
      filtered_count =+ 1
      return False
  else:
    vcf_count =+ 1
    return True
  return True

if __name__ == "__main__":
  main()
