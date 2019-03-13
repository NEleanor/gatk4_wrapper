"""Create the needed tool_data entry for the parse_gatk_json.py script"""

import sys

def main(file_name):
  json_file = open(file_name, 'r')
  #make a list for all the arguments
  arguments = []
  #initiate a tool name
  tool = file_name
  #iterate over the file line by line
  for line in json_file:
    #look for all lines that contain names of tools or arguments
    if "name" in line:
      #make sure we don't get the names of select options
      if "," in line:
        #look for the tool name which is the last one and the one without the --
        if "--" not in line:
          #split the line based on " and take the fourth piece 
          tool = line.split('"')[3]
        else:
          #if it is an argument then split the line as befor and add the argument name to the argument list
          argument = line.split('"')[3]
          arguments.append(argument)


  json_file.close()

  #a dictionary with all the macro names that go along with the arguments. macros must end in _cmd _inputs _pre or _outputs
  macro_dict = {"--help": ["all_cmd", "all_inputs", "log_outputs"], "--gatk-config-file": ["gatk_common_cmd", "gatk_common_inputs"], "--exclude-intervals": ["intervals_pre", "intervals_cmd", "intervals_inputs"], "--GA4GH_CLIENT_SECRETS": ["picard_cmd", "picard_inputs", "picard_outputs"], "--reference": ["ref_cmd", "ref_inputs"], "--REFERENCE_SEQUENCE": ["picard_ref_cmd", "ref_inputs"], "--read-filter": ["read_filter_cmd", "read_filter_inputs"], "--create-output-variant-md5": ["seq_dict_cmd", "seq_dict_inputs", "seq_dict_outputs"]}

  #initiate the strings with teh macros
  pre_commands = ""

  commands = ""

  inputs_params = ""

  output_params = ""

  #initate the list of macros
  macros = []

  #check each argument for corresponding macros
  for arg in arguments:
    if arg in macro_dict:
      macros = macros + macro_dict[arg]

  #checked each macro for its ending and add to the corresponding string
  for macro in macros:
    if "_cmd" in macro:
      if commands == "":
        commands += "'" + macro + "'"
      else:
        commands += ", '" + macro + "'"
    elif "_pre" in macro:
      if pre_commands == "":
        pre_commands += "'" + macro + "'"
      else:
        pre_commands += ", '" + macro + "'"
    elif "_inputs" in macro:
      if inputs_params == "":
        inputs_params += "'" + macro + "'"
      else:
        inputs_params += ", '" + macro + "'"
    elif "_outputs" in macro:
      if output_params == "":
        output_params += "'" + macro + "'"
      else:
        output_params += ", '" + macro + "'"

        
  #check if the tool name has an extra (Picard) and if so remove it
  if "(Picard)" in tool:
    tool = tool.split(" ")[0]
  
  #print the tool config and insert the needed tool name and macro strings
  print("'" + tool + "':\n\t{'output_fmt': {},\n\t'input_fmt':{},\n\t'pre_tmpls': [" + pre_commands + "],\n\t'post_tmpls': [" + commands + "'log_opts'],\n\t'pre_params':[" + inputs_params + "],\n\t'opt_params':[],\n\t'adv_params':[],\n\t'post_params': [],\n\t'output_params':[" + output_params + "]},")



if __name__ == "__main__":
    args = sys.argv
    main(args[1])