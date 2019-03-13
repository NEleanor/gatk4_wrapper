"""Remove the unnecessary cheeta codes"""

import sys

def main(file_name):
  #open the wrapper
  wrapper_file = open(file_name, 'r')
  #get the first line to find the name
  first_line = wrapper_file.readline()
  #make a default name in case somthing goes wrong
  new_file_name = first_line
  #find the name of the tool and create the name of the wrapper
  if '"' in first_line:
    new_file_name = first_line.split('"')[1] + "_post.xml"
  #print the new wrapper name
  print(new_file_name)
  #create and open the file
  new_wrapper = open(new_file_name, 'w')
  #write the first line to the new file
  new_wrapper.write(first_line)
  #the list of arguments that are in the macros and need to be removed from the cheetah section
  arguments = ["version", "showHidden", "help", "arguments_file", "VERBOSITY", "verbosity", "gatk_config_file", "gcs_max_retries", "use_jdk_deflater", "use_jdk_inflater", "interval_merging_rule", "interval_set_rule", "disable_read_filter", "disable_tool_default_read_filters", "read_filter", "create_output_variant_index", "add_output_sam_program_record", "add_output_vcf_command_line", "create_output_bam_index", "create_output_bam_md5", "create_output_variant_md5", "VALIDATION_STRINGENCY", "USE_JDK_INFLATER", "USE_JDK_DEFLATER", "TMP_DIR", "QUIET", "MAX_RECORDS_IN_RAM", "GA4GH_CLIENT_SECRETS", "CREATE_MD5_FILE", "CREATE_INDEX", "COMPRESSION_LEVEL", "REFERENCE_SEQUENCE", "OUTPUT", "SEQUENCE_DICTIONARY", "INPUT", "input", "reference", "output", "annotation", "annotation_group", "annotations_to_exclude", "intervals", "exclude_intervals", "read_index", "interval_padding", "interval_exclusion_padding", "output_prefix", "sequence_dictionary", "variant"]
  #initate the count to the high number
  count = 3
  #initate us outside the cheetah section
  in_cheetah = False
  #initiate the last line to an empty string
  last_line = ""
  #initiate keep_last to false so we don't prin the empty string
  keep_last = False
  #Loop over the wrapper file
  for line in wrapper_file:
    #initate keep to false
    keep = False
    # check if we are leaving the cheetah section
    if "<inputs>" in line:
      in_cheetah = False
    #check if we are in the cheetah section
    if in_cheetah:
      #loop over each argument in teh list of arguments to remove
      for argument in arguments:
        #find if the argument is in this line
        if "." + argument in line:
          #if the argument is in the line then set the count to remove this line and the next
          count = 1
      #if we are two away from a line with an argument then turn keep to true
      if count == 3:
        keep = True
      else:
        #if we aren't far enough away increase the count to one farther away
        count += 1
      if line.startswith("$"):
        #if the line starts with $ then it is a boolean and the line before and after should be removed so 
        #set the keep_last to false
        keep_last = False
        #set the count to 2 to get rid of the next one
        count = 2
    else:
      #if we are not in the cheetah section then keep the line
      keep = True
      #check if we are entering the cheetah section
      if "<command detect_errors" in line:
        in_cheetah = True
    #check if we are suppose to keep the last line and then write it if so
    if keep_last:
      new_wrapper.write(last_line)
    #check if we are in a boolean param
    if 'type="boolean"' in line:
      #check if we are in an improperly formatted default true boolean
      if 'checked="true"' in line and 'truevalue=""' not in line:
        #split the line by "
        line_split = line.split('"')
        #get the truevalue, add false to it and assing it to falsevalue
        new_false = line_split[7]
        line_split[7] = ""
        line_split[9] = new_false + " false"
        #initate the newline and put the split back together
        new_line = ""
        #initate a First check
        first = True
        for spot in line_split:
          # if first then skip the " and mark as no longer first
          if first:
            new_line = spot
            first = False
          else:
            #add the next spot and a "
            new_line += '"' + spot
          #assign the remade line to teh line
          line = new_line
    #remember the last line and it's keep status for the next interation. 
    last_line = line
    keep_last = keep
  #write the last line
  new_wrapper.write(last_line)
  #close the files
  wrapper_file.close()
  new_wrapper.close()
  
if __name__ == "__main__":
    args = sys.argv
    main(args[1])