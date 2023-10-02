import os
import re
import sys

# Define the regex pattern to match C++ variable declarations and their values
variable_pattern = r'\b(?:int|double|string)\s+(\w+)\s*=\s*([^;]+);'

# Initialize a flag to check if we are within the desired block
inside_block = False

# Initialize an empty dictionary to store variable names and values
variable_dict = {}

# Specify the input .cpp file path



def save_pram(file_path=''):
    '''
    python script to save cpp variables between\n
    //params-begin and //params-end as params.csv file\n
    Usage:
    give filename/path in the following priority
        1) call the function with filename as argument
            eg: save_pram('M6.cpp')
        2) run the main python file (as of now reg.py) followed by name 
            eg: python reg.py M6.cpp
        3) It will find the latest .cpp file in current working directory

    Example: \n
        cpp file: \n
            //params-begin \n   
            string child="cute"; \n
            int man=1;\n            
            double female=2.0;\n 
            //params-end\n
        params.csv: \n
            child,"cute"\n
            man,\n
            female,2.0\n
        '''
    if file_path:
        pass
    elif len(sys.argv)==2:
        file_path=sys.argv[1]
    else:
        files_in_directory = os.listdir()
        # Find the first ".cpp" file in the list
        cpp_files = [file for file in files_in_directory if file.endswith('.cpp')]
        if cpp_files[0]:
            cpp_files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
            file_path=cpp_files[0]
        else:
            print("No '.cpp' files found in the current directory.")
    print(file_path)
    # Read the file line by line and process it
    with open(file_path, 'r') as file:
        for line in file:
            # Check if we are entering the block
            if "//params-begin" in line:
                inside_block = True
                continue

            # Check if we are exiting the block
            if "//params-end" in line:
                inside_block = False
                break

            # If we are inside the block, try to match variable declarations and values
            if inside_block:
                match = re.match(variable_pattern, line.strip())
                if match:
                    variable_name = match.group(1)
                    variable_value = match.group(2)
                    variable_dict[variable_name] = variable_value

    # Print the extracted variable names and values
    with open("params.csv","w")as f:
        for variable_name, variable_value in variable_dict.items():
            f.write(f"{variable_name},{variable_value}\n")

if __name__=="__main__":
    save_pram()

