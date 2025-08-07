import re
import os
import math
import traceback
import subprocess
from itertools import combinations


atomsets = {}
#frist build a dictionary of all the atoms in a mmp using their index, proton_num, and xyz values also adding part number

#this dictionary holds the atoms from the mmp
atoms = {}
#this dictionary has the atom type data including atom type for output
output_atoms = {}
part_num = 0
atom_count = 0
atomset_name_list = []  #holds the names of the atomsets. If there are duplicate names, we'll need to append some suffix to one. 

Add_Topology = 1

#out_dir='C:/Users/Tom/Desktop/Blender Game Engine/assembly_test' # You can use forward slashes on Windows
out_dir="C:/LAMMPS/cantilever_positional_var/" # You can use forward slashes on Windows
data_file_out='temp' #this is a file that holds all the atoms before the header info is written to the final file which is another name down below
output_data=os.path.join(out_dir,data_file_out)
data_file_name = os.path.join(out_dir,"cantilever_positional_var.data")
mmp_file = "C:/LAMMPS/cantilever_positional_var/cantilever_positional_var.mmp"
input_file = open(mmp_file, 'r') 
igot = input_file.readlines()

for line_raw in igot:
    if line_raw.find("(Clipboard)") > -1:
        break    
    if line_raw.find("atomset") > -1 and "forward_ref" not in line_raw:
            #This section creates a dictionary for all the atomsets. It is possible to just use atomsets and not jigs. 
            jig_name_raw = re.findall(r'\(.*?\)', line_raw)
            atom_set_name = jig_name_raw[0].replace('(','').replace(')','')
            atomset_name_list.append(atom_set_name)
            
            if atom_set_name in atomset_name_list: #this allows multiple atomsets with the same name
                atom_set_name_X = atom_set_name + str(len(atomset_name_list))
            else:
                pass
            
            print("name list", atom_set_name)
            clean_line_1 = line_raw.replace('(0, 0, 0)', '')
            clean_line = clean_line_1.strip('\n')
            
            s1 = [int(i) for i in clean_line.split() if i.isdigit()] #this is the index of the atoms in the atomset           
            
            atomsets[atom_set_name_X] = s1
            

    line = line_raw.replace('(','').replace(')','')  
    
    line = line.replace(',','')
   
   #creates master dictionary of atoms from mmp
    if line.find("mol") > -1 and line.find("def") > -1:
        part_num = part_num + 1
    if line.find("atom ") > -1 and line.find("def") > -1:   

        mmp_atoms =  line.split()
        
        #print(mmp_atoms)
        
        atoms[mmp_atoms[1]] = {'part_num':part_num, 'proton_num':mmp_atoms[2], 'x_position':mmp_atoms[3], 'y_position':mmp_atoms[4], 'z_position':mmp_atoms[5] }   #this creates a nested dictionary for each atom
        
        atom_count = atom_count + 1
        
           
        

################################ATOM SETS############################################
#Change the atomset below if you have renamed them in the mmp
num_atom_types_and_mass = []                     #just use atoms sets and move away from jigs
try:
    for key, value in atomsets.items():        
        X_SUM = 0
        Y_SUM = 0
        Z_SUM = 0
        ATOM_NUM = 0
        
        for atoms_in_atomset in range(0,len(value)):            
            atomset_index = value[atoms_in_atomset]
            this_atoms_partnum = atoms[str(atomset_index)].get('part_num')
            this_atoms_proton_num = atoms[str(atomset_index)].get('proton_num')
            this_atoms_xpos = float(atoms[str(atomset_index)].get('x_position'))/1000
            this_atoms_ypos = float(atoms[str(atomset_index)].get('y_position'))/1000
            this_atoms_zpos = float(atoms[str(atomset_index)].get('z_position'))/1000   
            
            X_SUM = X_SUM + this_atoms_xpos
            Y_SUM = Y_SUM + this_atoms_ypos
            Z_SUM = Z_SUM + this_atoms_zpos
            ATOM_NUM = ATOM_NUM  +1 
            
            
            if "anchor" in key:
               
                
                if this_atoms_proton_num == "6":       
                        element = "C"
                        atom_type = "3"
                        atom_mass = 12.0109997
                        
                        if atom_type not in num_atom_types_and_mass:
                            num_atom_types_and_mass.append(atom_type)
                            num_atom_types_and_mass.append(str(atom_mass))
                            
                else:
                        element = "H"
                        atom_type = "4" 
                        atom_mass = 1.00796998
                        if atom_type not in num_atom_types_and_mass:
                            num_atom_types_and_mass.append(atom_type)
                            num_atom_types_and_mass.append(str(atom_mass))  
            
            elif "tip" in key:
                
                if this_atoms_proton_num == "6":       
                        element = "C"
                        atom_type = "5"
                        atom_mass = 12.0109997
                        if atom_type not in num_atom_types_and_mass:
                            num_atom_types_and_mass.append(atom_type)
                            num_atom_types_and_mass.append(str(atom_mass))                        
                else:
                        element = "H"
                        atom_type = "6"    
                        atom_mass = 1.00796998
                        if atom_type not in num_atom_types_and_mass:
                            num_atom_types_and_mass.append(atom_type)
                            num_atom_types_and_mass.append(str(atom_mass))                        
            elif key == "top":
                if this_atoms_proton_num == "6":       
                        element = "C"
                        atom_type = "7"
                        atom_mass = 12.0109997
                        if atom_type not in num_atom_types_and_mass:
                            num_atom_types_and_mass.append(atom_type)
                            num_atom_types_and_mass.append(str(atom_mass))                        
                else:
                        element = "H"
                        atom_type = "8"   
                        atom_mass = 1.00796998
                        if atom_type not in num_atom_types_and_mass:
                            num_atom_types_and_mass.append(atom_type)
                            num_atom_types_and_mass.append(str(atom_mass))
            elif key == "C60":
                if this_atoms_proton_num == "6":       
                        element = "C"
                        atom_type = "9"
                        atom_mass = 12.0109997
                        if atom_type not in num_atom_types_and_mass:
                            num_atom_types_and_mass.append(atom_type)
                            num_atom_types_and_mass.append(str(atom_mass))                        
                else:
                        element = "H"
                        atom_type = "7"   
                        atom_mass = 1.00796998
                        if atom_type not in num_atom_types_and_mass:
                            num_atom_types_and_mass.append(atom_type)
                            num_atom_types_and_mass.append(str(atom_mass))              
            elif key == "rotgroup":
                if this_atoms_proton_num == "6":       
                        element = "C"
                        atom_type = "6"
                        atom_mass = 12.0109997
                        if atom_type not in num_atom_types_and_mass:
                            num_atom_types_and_mass.append(atom_type)
                            num_atom_types_and_mass.append(str(atom_mass))                        
                else:
                        element = "H"
                        atom_type = "7"   
                        atom_mass = 1.00796998
                        if atom_type not in num_atom_types_and_mass:
                            num_atom_types_and_mass.append(atom_type)
                            num_atom_types_and_mass.append(str(atom_mass)) 
            
            #print( atomset_index, element, atom_type, this_atoms_xpos, this_atoms_ypos, this_atoms_zpos )
            #atom_line = str(atomset_index) + " 1 "  + " " +  " " + atom_type + " " + str(this_atoms_xpos) + " " +  str(this_atoms_ypos) + " " +  str(this_atoms_zpos) + "\n"
            #output_file.write(atom_line) 
            
            output_atoms[atomset_index] = {"mol_tag":this_atoms_partnum, "atom type":atom_type, "x cor":this_atoms_xpos, "y cor":this_atoms_ypos, "z cor":this_atoms_zpos}
            #remove it from the master atom dictionary
            atoms.pop(str(atomset_index))     
    
        x_centroid = round(X_SUM/ATOM_NUM, 3)
        y_centroid = round(Y_SUM/ATOM_NUM, 3)
        z_centroid = round(Z_SUM/ATOM_NUM, 3)
        
        print(key, "centroid", x_centroid, y_centroid, z_centroid)
    
    #num_atom_types = len(num_atom_types_and_mass)/2
    #print(num_atom_types_and_mass, "num of atom types", num_atom_types)    
    #handle the atoms that are left after removing the jig and atomset atoms 
except Exception:
    print(traceback.format_exc())
    pass
#############################THE REST OF THE ATOMS##############################

for atom_id in atoms.items():
    atom_index = int(atom_id[0])
    this_atoms_partnum = atoms[str(atom_id[0])].get('part_num')
    this_atoms_proton_num = atoms[str(atom_id[0])]['proton_num']
    this_atoms_xpos = float(atoms[str(atom_id[0])].get('x_position'))/1000
    this_atoms_ypos = float(atoms[str(atom_id[0])].get('y_position'))/1000
    this_atoms_zpos = float(atoms[str(atom_id[0])].get('z_position'))/1000  

    if this_atoms_proton_num == "6":       
            element = "C"
            atom_type = "1"
            atom_mass = 12.0109997
            if atom_type not in num_atom_types_and_mass:
                num_atom_types_and_mass.append(atom_type)
                num_atom_types_and_mass.append(str(atom_mass))                 
    else:
            element = "H"
            atom_type = "2"
            atom_mass = 1.00796998
            if atom_type not in num_atom_types_and_mass:
                num_atom_types_and_mass.append(atom_type)
                num_atom_types_and_mass.append(str(atom_mass))            
    
    output_atoms[atom_index] = {"mol_tag":this_atoms_partnum, "atom type":atom_type, "x cor":this_atoms_xpos, "y cor":this_atoms_ypos, "z cor":this_atoms_zpos}
    ##print(atom_index, element, atom_type, this_atoms_xpos, this_atoms_ypos, this_atoms_zpos )
    ##atom_line = str(atom_index) + " 1 "  + " " +  " " + atom_type + " " + str(this_atoms_xpos) + " " +  str(this_atoms_ypos) + " " +  str(this_atoms_zpos) + "\n"
    ##output_file.write(atom_line)  
num_atom_types = len(num_atom_types_and_mass)/2
print(num_atom_types_and_mass, "num of atom types", num_atom_types)   

#we've extracted differen atomset groups from the mmp e.g. atomsets and assigned a LAMMPS atomtype. Then we added those atoms to a new dictionary that gets sorted by index number. Sorting it probably not necessary, but we will do it anyway. 
lammps_atom_data_out = dict(sorted(output_atoms.items(), key=lambda item: item[0]))
for key, value in lammps_atom_data_out.items():
    #print(key, '->', value)
    #print(key)
    #print(lammps_atom_data_out[key].get('atom type'))
    
    atom_index_write = str(key)
    molecule_tag_write = str(lammps_atom_data_out[key].get('mol_tag')) #this is part number
    atom_type_write = str(lammps_atom_data_out[key].get('atom type'))
    x_cor_write = str(lammps_atom_data_out[key].get('x cor'))
    y_cor_write = str(lammps_atom_data_out[key].get('y cor'))
    z_cor_write = str(lammps_atom_data_out[key].get('z cor'))
                      
    atom_data_line = atom_index_write + "   " + molecule_tag_write + "   "  + atom_type_write + "   " + x_cor_write + "000"+  "   " + y_cor_write + "000"+ "  " + z_cor_write + "000" + "\n" #the "000" were added for a quick test, needs to be formated correctly
    with open(output_data,'a') as new_data:
        new_data.write(atom_data_line)

####get hi/lo xyz
Xs =  []
Ys = []
Zs = []
all_dim = []

written_data = open(output_data, 'r')

igot = written_data.readlines()

for line in igot:
    this = line.split()

    Xs.append(float(this[3]))
    Ys.append(float(this[4]))
    Zs.append(float(this[5]))
    

print("number of atoms in structure is : ", atom_count)
#print("lowest and highest Xs: ", min(Xs),max(Xs))
#print("lowest and highest Ys: ", min(Ys),max(Ys))
#print("lowest and highest Zs: ", min(Zs),max(Zs))

padding = 5

X_low =  min(Xs)
X_high = max(Xs)
all_dim.append(X_high)

Y_low =  min(Ys)
Y_high = max(Ys)
all_dim.append(Y_high)

Z_low =  min(Zs)
Z_high = max(Zs)
all_dim.append(Z_high)

#print("highest value", max(all_dim))
highest_value = max(all_dim)



if X_low < 0:
    X_low_pad = math.floor(X_low - padding)
else:
    X_low_pad = math.floor(X_low + padding)
    
if X_high < 0:
    X_high_pad = math.floor(X_high - padding)
else:
    X_high_pad = math.floor(X_high + padding)

if Y_low < 0:
    Y_low_pad = math.floor(Y_low - padding)
else:
    Y_low_pad = math.floor(Y_low + padding)
    
if Y_high < 0:
    Y_high_pad = math.floor(Y_high - padding)
else:
    Y_high_pad = math.floor(Y_high + padding)    


if Z_low < 0:
    Z_low_pad = math.floor(Z_low - padding)
else:
    Z_low_pad = math.floor(Z_low + padding)
    
if Z_high < 0:
    Z_high_pad = math.floor(Z_high - padding)
else:
    Z_high_pad = math.floor(Z_high + padding)    


#print("XLOW_padding", X_low_pad, "XHIGH_padding", X_high_pad )     
written_data.close()
#####This was added after writing out the Atoms sections#####

comment = data_file_name + "\n"
num_atoms = " " + str(atom_count) + " atoms" + "\n"
atom_type_header = " " + str(int(num_atom_types)) + " atom types" + "\n"
xhilo = "   "+ str(X_low_pad) + "          "+ str(X_high_pad) + "      "+ "xlo xhi"+  "\n"
yhilo = "   "+ str(Y_low_pad) + "          "+ str(Y_high_pad) + "      "+ "ylo yhi"+  "\n"
zhilo = "   "+ str(Z_low_pad) + "          "+ str(Z_high_pad) + "      "+ "zlo zhi"+  "\n"
masses_header = " Masses" + "\n"
atoms_header = " Atoms " + "#molecular"+ "\n"


f = open(data_file_name, "a") 
#f = open(output_data, "a")
f.write(comment)
f.write("\n")
f.write(num_atoms)
f.write("\n")
f.write(atom_type_header)
f.write("\n")
f.write(xhilo)
f.write(yhilo)
f.write(zhilo)
f.write("\n")
f.write(masses_header)
f.write("\n")

for items in range(0, len(num_atom_types_and_mass), 2):
    f.write(" " + num_atom_types_and_mass[items] + " " +  num_atom_types_and_mass[items+1] + "\n")
f.write("\n")
f.write(atoms_header)
f.write("\n")


with open(output_data,'r') as atom_data:
    igot = atom_data.readlines()
    for line_raw in igot:
        f.write(line_raw)
f.close()

os.remove(output_data)
#################################################################################
if Add_Topology == 1:  
    print("adding topology")
    subprocess.run(["python", "add_lammps_bonds_7.py", out_dir, mmp_file])
    
    
    #out_dir='C:/Users/Tom/Desktop/Blender Game Engine/assembly_test' # You can use forward slashes on Windows
    out_dir="C:/LAMMPS/cantilever_positional_var/" # You can use forward slashes on Windows
    data_file_out='temp' #this is a file that holds all the atoms before the header info is written to the final file which is another name down below
    output_data=os.path.join(out_dir,data_file_out)
    data_file_name = os.path.join(out_dir,"cantilever_positional_var.data")
    mmp_file = "C:/LAMMPS/cantilever_positional_var/cantilever_positional_var.mmp"
    input_file = open(mmp_file, 'r') 
    
    
    #get the number of bonds
    bonds = open(out_dir+"Bonds.txt", 'r')
    igot_b = bonds.readlines()
    number_bonds = igot_b[-1].split()[0]
    print("number of bonds", number_bonds)
    bonds.close()
    #get number of angles
    angles = open(out_dir+"Angles_numbered.txt", 'r')
    igot_a = angles.readlines()
    number_angles = igot_a[-1].split()[0]
    print("number of angles", number_angles)
    angles.close()
    
    num_bonds_line = 3
    num_angles_line = 4
    bond_coef_line = 7
    angle_coef_line = 8
   
    file_path = "C:/LAMMPS/cantilever_positional_var/cantilever_positional_var.data"
    bond_path = "C:/LAMMPS/cantilever_positional_var/Bonds.txt"
    angle_path = "C:/LAMMPS/cantilever_positional_var/Angles_numbered.txt"

    
    with open(file_path, 'r') as file:
        lines = file.readlines()
        
    #insert the number of bonds
    bond_number_str = str(number_bonds)+" "+"bonds" + "\n"
    lines.insert(num_bonds_line, bond_number_str )
    #insert number of angles with carrage return after
    angle_number_str = str(number_angles)+" "+"angles" +"\n"
    lines.insert(num_angles_line, angle_number_str )
    
    #insert number of bond and angle coeffients 
    bond_coef_str =" "+ "2 bond types" + "\n" # needs to be changed it there are no hydrogen, eventually it should count the number in the bonds list
    lines.insert(bond_coef_line, bond_coef_str )
   
    angle_coef_str =" "+ "3 angle types"  +"\n"  # needs to be changed it there are no hydrogen, eventually it should count the number in the angles list
    lines.insert(angle_coef_line , angle_coef_str )    
        
    with open(file_path, 'w') as file:
        file.writelines(lines)   
        
    with open(file_path, 'a') as file:
        blank_line = '\n'
        
        file.writelines(blank_line)
        file.writelines(" Bonds")   
        file.writelines(blank_line)
        file.writelines(blank_line)
        with open(bond_path, 'r') as bond_file:
            bonds = bond_file.readlines()
            for bond in bonds:
                file.writelines(bond)
      
        #angles section
        file.writelines(blank_line)
        file.writelines(blank_line)       
        
        file.writelines(blank_line)
        file.writelines(" Angles")   
        file.writelines(blank_line)
        file.writelines(blank_line)
        with open(angle_path, 'r') as angles_file:
            angles= angles_file.readlines()
            for angle in angles:
                file.writelines(angle)    

print("LAMMPS DATA FILE COMPLETED. CHECK THE BOX SIZE")    
    
 
        
  
    
    

        
    
    

    
    
    
    


