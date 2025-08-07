import os
import sys
from itertools import combinations

##################################
atoms_and_bonds = {} 

if len(sys.argv) > 1:
            project_path = sys.argv[1]
            mmp_file =  sys.argv[2]            
else:

            project_path = "C:\\LAMMPS\\cantilever_positional_var\\SAMSON topo\\"  
            mmp_file = project_path + "cantilever_positional_var.mmp"

#CHS single carbon hydrogen bond, bond number 2
#CCS single carbon carbon bond, bond number 1

with open(mmp_file , 'r') as file:
            igot = file.readlines()
            for atom_line_num, line in enumerate(igot, 1):
                        
                        if line.find("atom") > -1:
                                    atom_line_info = line.split()
                                    atom_index = atom_line_info[1]
                                    atom_element =  atom_line_info[2]
                                    #print("atom info line", line)
                                    if atom_element == ("(6)"): 
                                                atoms_and_bonds[atom_index] = {'element':6, 'valence':4, 'index_bond1':'0', 'type_bond1':'0', 'index_bond2':'0', 'type_bond2':'0', 'index_bond3':'0', 'type_bond3':'0', 'index_bond4':'0', 'type_bond4':'0' } #this needs to be index_bonded_to, bondtype,
                                    else:
                                                atoms_and_bonds[atom_index] = {'element':1, 'valence':1, 'index_bond1':'0', 'type_bond1':'0'} 

                        if line.find("bond") > -1:
                                    #print(line)
                                    bonded_to_list = []
                                    bond_info = line.strip("bond").split()
                                    #print("bond info", bond_info)
                                    bond_type = bond_info[0]
                                    #print("bond type", bond_type)
                                    bonded_to_list.append(bond_info[1])  
                                    #print("bonded to list first time", bonded_to_list)
                                    for otherbond in range(2, 5):
                                                
                                                try:
                                                            bonded_to_list.append(bond_info[otherbond])
                                                            #print("bonded to list loop", bonded_to_list)
                                                except:                                                            
                                                            break
                                  
                                    for other_bonds in range(len(bonded_to_list) + 1):                                    
                                                #print("atom index ", atom_index," element ", atom_element)
                                                
                                                try:
                                                
                                                            bonded_to = bonded_to_list[other_bonds]
                                                except:
                                                            bonded_to = None
                                                            break
                                                
                                                #print("bond type:", bond_type, "bonded to atom index1: ", bonded_to)
                                                                                              
                                                key_list = list(atoms_and_bonds.keys())
                                                
                                                bonded_to_index = int(bonded_to)-1   #the list starts at zero, just a double check anyway that it is opening the write dictionary                                    

                                                #print("atom ", key_list[bonded_to_index], " is the dictionary to add the bonded to info to", atoms_and_bonds[bonded_to])
                                                
                                                ##for every atom index, open the dictionary that that atom is bonded to and add each bond type and the bonded to information, 
                                                                                               
                                                num_bonds = atoms_and_bonds[bonded_to]['valence']
                                                
                                                for valences in range(1,num_bonds+1):
                                                            bond_to_update = "index_bond"+str(valences)
                                                            #print(atoms_and_bonds[bonded_to].get(bond_to_update))
                                                                                                            
                                                            insert = atoms_and_bonds[bonded_to].get(bond_to_update)
                                                            if insert == '0':
                                                                        #print("bond_to_update", bond_to_update)
                                                                        #print("atom index", atom_index)
                                                                        atoms_and_bonds[str(bonded_to)][bond_to_update] = atom_index  
                                                                        
                                                                        #print("atom index", atom_index)
                                                                        break
                                                            elif insert !=0:
                                                                        continue
 
                                                for valences in range(1,num_bonds+1):          
                                                            type_to_update = "type_bond"+str(valences)  
                                                            insert_type = atoms_and_bonds[bonded_to].get(type_to_update)
                                                            if insert_type == '0':
                                                                        #print(atoms_and_bonds[bonded_to]['element'], "below", atom_element, atom_index)
                                                                        if atom_element == "(6)": #and bond_type == "1":   
                                                                                    if atoms_and_bonds[bonded_to]['element'] == 6:
                                                                                                atoms_and_bonds[str(bonded_to)][type_to_update] = "1" # "CCS"  #carbon-carbon single bond
                                                                                    if atoms_and_bonds[bonded_to]['element'] == 1:
                                                                                                atoms_and_bonds[str(bonded_to)][type_to_update] = "2" #carbon-hydrogen single bond
                                                                        else:
                                                                                    
                                                                                    atoms_and_bonds[str(bonded_to)][type_to_update] = "2"  # "CHS"  #carbon-hydrogen single bond
                                                                        break
                                                            elif type_to_update != 0:
                                                                        continue

                                  
#####un-comment out to print the connection dictionary#################################
#print('PRINT CONNECTION DICTIONARY')
#for atom_index in atoms_and_bonds:            
            #print('atom_index = ',atom_index)
            #for inner_key in atoms_and_bonds[atom_index]:
                        #print((inner_key,atoms_and_bonds[atom_index][inner_key]))                                      
#####################################################################################                        

##write out bonds. LAMMPS format is Bond number, bond coeff type, first index, last atom index e.g. 1 5 1 2
bond_number = 0
for outer_key in atoms_and_bonds:
            me_element = atoms_and_bonds[outer_key]['element']
            num_bonds = atoms_and_bonds[outer_key]['valence']
            
            for valences in range(1,num_bonds+1):
                        last_atom_index = "index_bond"+str(valences)
                        bond_coeff_type = "type_bond" + str(valences)
                        if atoms_and_bonds[outer_key][bond_coeff_type] != "0":
                                    bond_number =  bond_number + 1
                        
                                    #print("bond number", bond_number,  "bond coefficient", atoms_and_bonds[outer_key][bond_coeff_type], " frist index", outer_key, "last atom index", atoms_and_bonds[outer_key][last_atom_index ])
                                    bond_line_to_write = str(bond_number) + " " + str(atoms_and_bonds[outer_key][bond_coeff_type]) + " " + str(outer_key) + " " + str(atoms_and_bonds[outer_key][last_atom_index ]) + "\n"
                                    with open(project_path + 'Bonds.txt', 'a') as file:
                                                file.write(bond_line_to_write)

####Angles
#LAMMPS DATA FORMAT 1 4 1 2 3: angle number, angle coreffecient number, atoms comprising the angle. LAMMPS docs on formating data confirms the middle atom is the corner i.e. KEY atom in the dictionary                                    
#H-CS-H, angle type 1
#CS-CS-H angle type 2
#CS-CS-CS angle type 3
#CHS single carbon hydrogen bond, bond number 2
#CCS single carbon carbon bond, bond number 1

def get_0_angles(end_carbon, outer_key, other_end_index, other_end_bondtype ):
            #deduce the angle corefficient from the key element and two corner bond type
                                                                       
            if  other_end_bondtype == "2" :
                        angle_type = "2" #"CCH"
                        #print("CCH"," angle type ", angle_type)  
            else:
                        angle_type = "3" #"CCC"
                        #print("CCC"," angle type ", angle_type)
            angle_line_to_write = angle_type + " " + end_carbon + " " + outer_key + " " + other_end_index + "\n"
            #print(angle_line_to_write)
            
            with open(project_path + 'Angles.txt', 'a') as file:
                        file.write(angle_line_to_write)                                                                          
            

for outer_key in atoms_and_bonds:            
            key_element = atoms_and_bonds[outer_key]['element']
            num_bonds = atoms_and_bonds[outer_key]['valence'] 
            #print(outer_key, key_element, num_bonds)
            zero_count = False
            if key_element != 1: #not hydrogen atoms. #C5 has a 0 in its dic record for its bond with C1. This is causing angles 651, 751, 851 to be skipped. If a C has a 0 in its records, need to find the other carbon that has its index
                        element_letter = "C"
                        #print("middle atom is index",outer_key, "and the element is", element_letter)
                        for outer_valences in range(1,5): 
                                    #print("atom index", outer_key, "outer_valence", outer_valences)
                                    corner_1_index = atoms_and_bonds[outer_key]['index_bond' + str(outer_valences)]
                                    corner_1_bondtype = atoms_and_bonds[outer_key]['type_bond' + str(outer_valences)]
                                    if corner_1_index == "0" and corner_1_bondtype == "0":
                                                
                                                #print("big zero on carbon with atom index:", outer_key )
                                                #print("let's find the other carbon atom its connected to.....")
                                                if zero_count == False:
                                                            zero_count = True
                                                            for atoms in atoms_and_bonds: 
                                                                        for bonded_atoms in range(1,5):
                                                                                    bond_to_search = "index_bond"+str(bonded_atoms)                                                                        
                                                                                    search = (atoms_and_bonds[atoms].get(bond_to_search))
                                                                                    if search == outer_key:
                                                                                                #print("match found in atoms_and_bonds dictionary", bond_to_search )
                                                                                                #print("these are my connect records since I am atom index", outer_key, atoms_and_bonds[outer_key])
                                                                                                #print("I'm connected to atom with atom index", atoms, atoms_and_bonds[atoms])
                                                                                                #print("need to form angles with atoms in my connect record and atom", atoms)
                                                                                                middle_atom = outer_key
                                                                                                end_carbon = atoms
                                                                                                for connected_atoms in range(1,5):
                                                                                                            other_end_index = atoms_and_bonds[outer_key]['index_bond' + str(connected_atoms)]
                                                                                                            other_end_bondtype = atoms_and_bonds[outer_key]['type_bond' + str(connected_atoms)]
                                                                                                            if other_end_index != "0":
                                                                                                                        #print(" ")
                                                                                                                        #print(other_end_index, other_end_bondtype)
                                                                                                                        #one end and middle atom will always be C, so only need their index number not bondtype
                                                                                                                        #the second end atom might be C or H?, so test that with the bondtype to assign the angletype
                                                                                                                        get_0_angles(end_carbon, outer_key, other_end_index, other_end_bondtype )
                                                                                   
                                                            


                                                
                                    else:
                                                pass
                                                #print("                    no 0s")
                                    
                                    #print("j", outer_valences, "corner one atom index", corner_1_index, "corner 1 bond type", corner_1_bondtype)
                                    for  inner_valences in range(outer_valences+1,num_bonds+1):
                                                #print("k", inner_valences)
                                                corner_2_index = atoms_and_bonds[outer_key]['index_bond' + str(inner_valences)]
                                                corner_2_bondtype = atoms_and_bonds[outer_key]['type_bond' + str(inner_valences)] 
                                                #print("k", inner_valences, "corner two atom index", corner_2_index, "corner 2 bond type", corner_2_bondtype)
                                                
                                                #deduce the angle corefficient from the key element and two corner bond type
                                                if key_element == 6 and corner_1_bondtype == "2"  and corner_2_bondtype == "2":
                                                            angle_type = "1"  #"HCH"
                                                            #print("HCH", " angle type ", angle_type)                                                
                                                if key_element == 6 and corner_1_bondtype == "1"  and corner_2_bondtype == "2" or corner_1_bondtype == "2"  and corner_2_bondtype == "1" :
                                                            angle_type = "2" #"CCH"
                                                            #print("CCH"," angle type ", angle_type)  
                                                if key_element == 6 and corner_1_bondtype == "1"  and corner_2_bondtype == "1":
                                                            angle_type = "3" #"CCC"
                                                            #print("CCC"," angle type ", angle_type)                                                    
                                                
                                                
                                               
                                                
                                                
                                                #create the angle string to write to the datafile
                                                #some test to see if a corner index is zero and if so pass on through without writting the angle
                                                
                                                if  corner_2_index != "0":
                                                            angle_line_to_write = angle_type + " " + corner_1_index + " " + outer_key + " " + corner_2_index + "\n"
                                                            #print(angle_line_to_write)
                                                            
                                                            with open(project_path + 'Angles.txt', 'a') as file:
                                                                        file.write(angle_line_to_write)                                                              
                                                else:
                                                            pass
              

def get_more0s_angles(outer_key, zero_count, zero_list):
            atom_angle_list = [ ]
            
            for index in range(0, len(zero_list), 2):
                        atom_angle_list.append(zero_list[index])
            print(outer_key, zero_count, zero_list)
            if zero_count == 2:
                        print(outer_key, zero_count, zero_list)
                        first_zero = zero_list[0]
                        first_bond = zero_list[1]
                        second_zero = zero_list[2]
                        second_bond = zero_list[3]
                       
                        if  first_bond == "2" and second_bond == "2" :
                                    angle_type = "1" #"HCH"
                                    #print("CCH"," angle type ", angle_type)  
                        elif first_bond == "1" and second_bond == "2":
                                    angle_type = "2" #"CCH"
                                    #print("CCC"," angle type ", angle_type)
                        else:
                                    angle_type = "3" #"CCC"   
                                    
                        angle_line_to_write = angle_type + " " + first_zero + " " + outer_key + " " + second_zero + "\n"
                        with open(project_path + 'Angles.txt', 'a') as file:
                                    print("writing a 2 zero")
                                    file.write(angle_line_to_write)                         
                        
            
            
            if zero_count == 3:
                        first_bond = 0
                        second_bond = 0
                        #print(outer_key, zero_count, zero_list)
                        couplets = list(combinations(atom_angle_list, 2))                                               
                        for couplet in couplets:
                                    #print(couplet)
                                    for pairs in range(2):
                                                pair_index = zero_list.index(couplet[pairs])
                                                #print(pair_index)
                                                bond_index = int(pair_index) + 1
                                                if first_bond == 0:
                                                            #print(bond_index)
                                                            first_bond = zero_list[bond_index]
                                                else:
                                                            second_bond = zero_list[bond_index]
                                    
                                    first_zero = couplet[0]
                                    second_zero = couplet[1]
                                    
                                    #print(first_zero, first_bond, outer_key, second_bond, second_zero )
                                    
                                    if  first_bond == "2" and second_bond == "2" :
                                                angle_type = "1" #"HCH"
                                                #print("CCH"," angle type ", angle_type)  
                                    elif first_bond == "1" and second_bond == "2":
                                                angle_type = "2" #"CCH"
                                                #print("CCC"," angle type ", angle_type)
                                    else:
                                                angle_type = "3" #"CCC"    
                                    angle_line_to_write = angle_type + " " + first_zero + " " + outer_key + " " + second_zero + "\n"
                                    with open(project_path + 'Angles.txt', 'a') as file:
                                                print("writing a 3 zero")
                                                file.write(angle_line_to_write)                                 
                        
            if zero_count == 4:
                        first_bond = 0
                        second_bond = 0
                        print(outer_key, zero_count, zero_list)
                        couplets = list(combinations(atom_angle_list, 2))                                               
                        for couplet in couplets:
                                    print(couplet)
                                    for pairs in range(2):
                                                pair_index = zero_list.index(couplet[pairs])
                                                #print(pair_index)
                                                bond_index = int(pair_index) + 1
                                                if first_bond == 0:
                                                            ##print(bond_index)
                                                            first_bond = zero_list[bond_index]
                                                else:
                                                            second_bond = zero_list[bond_index]
                                    
                                    first_zero = couplet[0]
                                    second_zero = couplet[1]
                                    
                                    #print(first_zero, first_bond, outer_key, second_bond, second_zero )
                                    
                                    if  first_bond == "2" and second_bond == "2" :
                                                angle_type = "1" #"HCH"
                                                print("HCH"," angle type ", angle_type)  
                                    elif first_bond == "1" and second_bond == "2":
                                                angle_type = "2" #"CCH"
                                                print("CCH"," angle type ", angle_type)
                                    else:
                                                angle_type = "3" #"CCC"    
                                    angle_line_to_write = angle_type + " " + first_zero + " " + outer_key + " " + second_zero + "\n"
                                    print(angle_line_to_write)                        
            
            
                                    with open(project_path + 'Angles.txt', 'a') as file:
                                                print("writing a 4 zero")
                                                file.write(angle_line_to_write)   



#handle carbon atoms that have two zeros (or more?) in their connect records on their own so as not to over complicate the above. There are more, carbons with 3 and even 4 zeros. 
#lets try and change this to use a list that hold 2,3,4 zeros and write the angle combos for everything in the list

#This will handle carbons with more than two zero. 
for outer_key in atoms_and_bonds: 
            zero_count = 0
            #print("initializing zero list")
            zero_list = [ ] 
            key_element = atoms_and_bonds[outer_key]['element']
            num_bonds = atoms_and_bonds[outer_key]['valence']
            if key_element != 1: #not hydrogen atoms. #C5 has a 0 in its dic record for its bond with C1. This is causing angles 651, 751, 851 to be skipped. If a C has a 0 in its records, need to find the other carbon that has its index
                        element_letter = "C"
                        #print("middle atom is index",outer_key, "and the element is", element_letter)
                        for outer_valences in range(1,5): 
                                    #print("atom index", outer_key, "outer_valence", outer_valences)
                                    corner_1_index = atoms_and_bonds[outer_key]['index_bond' + str(outer_valences)]
                                    corner_1_bondtype = atoms_and_bonds[outer_key]['type_bond' + str(outer_valences)]
                                    if corner_1_index == "0" and corner_1_bondtype == "0":
                                                zero_count = zero_count + 1
                                                #print("big zero on carbon with atom index:", outer_key) 
                        if zero_count > 1: #atoms with only one zero are handled in the above code                                   
                                    #print(outer_key, "num zeros", zero_count) 
                                    
                                    for atoms in atoms_and_bonds: 
                                                for bonded_atoms in range(1,5):
                                                            bond_to_search = "index_bond"+str(bonded_atoms)                                                                        
                                                            search = atoms_and_bonds[atoms].get(bond_to_search)                                                                       
                                                            if search == outer_key:   
                                                                        #print("match found in atoms_and_bonds dictionary", bond_to_search )
                                                                        bond_to_get = "type_bond" +str(bonded_atoms)
                                                                        bond_type = atoms_and_bonds[atoms].get(bond_to_get)
                                                                        zero_list.append(atoms)
                                                                        zero_list.append(bond_type)
                                    #print(outer_key, "num zeros", zero_count, zero_list)
                                    get_more0s_angles(outer_key, zero_count, zero_list)

                                     
      

with open(project_path + 'Angles.txt', 'r') as file:
            angle_number  =  1
            igot = file.readlines()    
            for line in igot:    
                        angle_info = line.split()
                        #print(angle_info)
                        with open(project_path + 'Angles_numbered.txt', 'a') as file2:
                                    write_angle = str(angle_number) + " " + angle_info[0] + " " + angle_info[1] + " " + angle_info[2] + " " + angle_info[3] + "\n"
                                    file2.write(write_angle)
                                    angle_number  =  angle_number + 1

os.remove(project_path + 'Angles.txt')


#Note for if I ever want to add a dihedral section. For each atom (outter key) go through what atoms it is bonded to (inner key),
# for each of those atoms, use their index number (outter_key-2) to look at their atoms_and_bonds dictionary to get what next over atom they connect to
#use this to identify the dihedral angles



