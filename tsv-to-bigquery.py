import os
import re
import sys

def main():
	  """Entry point to the script."""

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
      accumulate_block(fields)

    line = file_handle.readline()

  # Emit the final block, if applicable
  emit_block(sample_id)


if __name__ == "__main__":
	main()
	print "Done!"	