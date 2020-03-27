#Modified version 2/28/20

'''
	Up to 88 rxns per single run and 24 rxns for triplicate of each combination of DNA assembly
	Multichannel p300, single channel p10
    This protocol is optimized for maximum walkaway time and percision
	'''

import time
import math

from opentrons import robot, instruments, labware, modules

'''
COLD_BLOCK = '96-PCR-tall-cold-block'
try:
	labware.create(
				   COLD_BLOCK,
				   grid=(12, 8),
				   spacing=(9, 9),
				   diameter=5,
				   depth=15.4,
				   volume=200)
except:
	print("Using existing labware definition for {0}".format(COLD_BLOCK))
'''

num_rxns = len(combinations_to_make)
#num_plates = math.ceil(num_rxns/24) # Amount of agar plates need for plating the transformed cells

'''
	For this protocol, the biorad_96_wellplate_200ul_pcr, is placed on the top of the TempDeck alone without the Opentrons 96 well aluminum block.
	'''

# Load in Bio-Rad 96 Well Plate on temp deck for moclos, transformation, and outgrowth.
temp_deck = modules.load('tempdeck', '10')
reaction_plate = labware.load('biorad_96_wellplate_200ul_pcr', '10', share=True)
temp_deck.set_temperature(10)

# Load in 1 10ul tiprack and 2 300ul tipracks
tr_10 = [labware.load('opentrons_96_tiprack_10ul', '3')]#, labware.load('opentrons_96_tiprack_10ul', '6')]
tr_300 = [labware.load('opentrons_96_tiprack_300ul', '6'), labware.load('opentrons_96_tiprack_300ul', '9')]
#for i in range(0, 1):
#   tr_300.append(labware.load('tipone_96_tiprack_200ul', '9'))

# Load in pipettes
p10_single = instruments.P10_Single(mount='right', tip_racks=tr_10)
p300_multi = instruments.P300_Multi(mount='left', tip_racks=tr_300)

''' Need to provide the instructions for loading reagent'''
reagents_plate = labware.load('biorad_96_wellplate_200ul_pcr', '4', 'Reagents Plate')
ligase = reagents_plate.wells('H12') #MoClo
restriction_enzyme = reagents_plate.wells('G12') #MoClo
buffer = reagents_plate.wells('F12') # MoClo

'''
	This deck slot location is dedicated for the reaction plate after MoClo protocol is completed, so at the beginning of the protocol there isn't an actual plate existing in this slot location.
	'''
post_moclo_reaction_plate = labware.load('biorad_96_wellplate_200ul_pcr', '7', 'Post-MoClo Reaction Plate')

# Load in water, SOC, and wash trough (USA Scientific 12 Well Reservoir 22ml)
trough = labware.load('usascientific_12_reservoir_22ml', '5', 'Reagents trough')
water = trough.wells(0) #Well 1
wash_0 = trough.wells(1) #Well 2
wash_1 = trough.wells(2) #Well 3
soc = trough.wells(3) #Well 4
liquid_waste = trough.wells(4) #Well 5

# Load in up to 2 DNA plates (Bio-Rad 96 Well Plate 200ul PCR)
#plate_name = dna_plate_map_dict.keys() # because there should be only 1 input plate
#input_dna_plate = labware.load('biorad_96_wellplate_200ul_pcr', '1', 'Input DNA Plate')
dna_plate_dict = {}
for plate_name in dna_plate_map_dict.keys():
	dna_plate_dict[plate_name] = labware.load('biorad_96_wellplate_200ul_pcr', '1', 'Input DNA Plate')

# Load in 1 agar plate, same antibiotic for all plasmids is assumed (e-gelgol)
# Be careful labware e-gelgol is an old container type, and could eventually be removed by Opentrons!
agar_plate = labware.load('e-gelgol', '2', 'Agar Plate')

#available_deck_slots = ['11', '8']

#____________________________________Start the MoClo protocol_______________________________

'''
	For this protocol, the biorad_96_wellplate_200ul_pcr, is placed on the top of the TempDeck alone without the metal part.
	'''
# Add water, buffer, restriction enzyme, ligase, and buffer to 2x master mix (2xMM).
#Prepare 2xMM for 2-part assembly at 20uL total volume
for i in range(2):
	p10_single.pick_up_tip()
	p10_single.transfer((150), water.bottom(), reaction_plate.wells(95-i).bottom(0.5), new_tip='never') #Water
	p10_single.blow_out()
	p10_single.drop_tip()
	p10_single.pick_up_tip()
	p10_single.transfer(24, buffer.bottom(), reaction_plate.wells(95-i).bottom(0.5), new_tip='never') #Buffer (1uL/rxn)
	p10_single.mix(2, 10, reaction_plate.wells(95-i).bottom(0.5))
	p10_single.blow_out()
	p10_single.drop_tip()
	p10_single.pick_up_tip()
	p10_single.transfer(6, ligase.bottom(), reaction_plate.wells(95-i).bottom(0.5), new_tip='never') #Ligase (0.25uL/rxn)
	p10_single.mix(2, 10, reaction_plate.wells(95-i).bottom(0.5))
	p10_single.blow_out()
	p10_single.drop_tip()
	p10_single.pick_up_tip()
	p10_single.transfer(12, restriction_enzyme.bottom(), reaction_plate.wells(95-i).bottom(0.5), new_tip='never') #Restriction Enzyme  (0.5uL/rxn)
	p10_single.mix(2, 10, reaction_plate.wells(95-i).bottom(0.5))
	p10_single.blow_out()
	p10_single.drop_tip()

#Prepare 2xMM for 5-part assembly at 20uL total volume
for i in range(2):
	p10_single.pick_up_tip()
	p10_single.transfer((78), water.bottom(), reaction_plate.wells(93-i).bottom(0.5), new_tip='never') #Water
	p10_single.blow_out()
	p10_single.drop_tip()
	p10_single.pick_up_tip()
	p10_single.transfer(24, buffer.bottom(), reaction_plate.wells(93-i).bottom(0.5), new_tip='never') #Buffer (1uL/rxn)
	p10_single.mix(2, 10, reaction_plate.wells(93-i).bottom(0.5))
	p10_single.blow_out()
	p10_single.drop_tip()
	p10_single.pick_up_tip()
	p10_single.transfer(6, ligase.bottom(), reaction_plate.wells(93-i).bottom(0.5), new_tip='never') #Ligase (0.25uL/rxn)
	p10_single.mix(2, 10, reaction_plate.wells(93-i).bottom(0.5))
	p10_single.blow_out()
	p10_single.drop_tip()
	p10_single.pick_up_tip()
	p10_single.transfer(12, restriction_enzyme.bottom(), reaction_plate.wells(93-i).bottom(0.5), new_tip='never') #Restriction Enzyme  (0.5uL/rxn)
	p10_single.mix(2, 10, reaction_plate.wells(93-i).bottom(0.5))
	p10_single.blow_out()
	p10_single.drop_tip()

#Prepare 2xMM for 8-part assembly at 20uL total volume
for i in range(2):
	p10_single.pick_up_tip()
	p10_single.transfer((6), water.bottom(), reaction_plate.wells(91-i).bottom(0.5), new_tip='never') #Water
	p10_single.blow_out()
	p10_single.drop_tip()
	p10_single.pick_up_tip()
	p10_single.transfer(24, buffer.bottom(), reaction_plate.wells(91-i).bottom(0.5), new_tip='never') #Buffer (1uL/rxn)
	p10_single.mix(2, 10, reaction_plate.wells(91-i).bottom(0.5))
	p10_single.blow_out()
	p10_single.drop_tip()
	p10_single.pick_up_tip()
	p10_single.transfer(6, ligase.bottom(), reaction_plate.wells(91-i).bottom(0.5), new_tip='never') #Ligase (0.25uL/rxn)
	p10_single.mix(2, 10, reaction_plate.wells(91-i).bottom(0.5))
	p10_single.blow_out()
	p10_single.drop_tip()
	p10_single.pick_up_tip()
	p10_single.transfer(12, restriction_enzyme.bottom(), reaction_plate.wells(91-i).bottom(0.5), new_tip='never') #Restriction Enzyme  (0.5uL/rxn)
	p10_single.mix(2, 10, reaction_plate.wells(91-i).bottom(0.5))
	p10_single.blow_out()
	p10_single.drop_tip()

# Add master mix to each wells that would contain rxn
# Add 2xMM to 2-part assembly
p10_single.pick_up_tip()
for i in range(12):
	p10_single.transfer(16, reaction_plate.wells(95).bottom(), reaction_plate.wells(i).bottom(0.3), new_tip='never')
	p10_single.blow_out()
p10_single.drop_tip()

p10_single.pick_up_tip()
for i in range(12):
	p10_single.transfer(16, reaction_plate.wells(94).bottom(), reaction_plate.wells(12+i).bottom(0.3), new_tip='never')
	p10_single.blow_out()
p10_single.drop_tip()

# Add 2xMM to 5-part assembly
p10_single.pick_up_tip()
for i in range(12):
	p10_single.transfer(10, reaction_plate.wells(93).bottom(), reaction_plate.wells(24+i).bottom(0.3), new_tip='never')
	p10_single.blow_out()
p10_single.drop_tip()

p10_single.pick_up_tip()
for i in range(12):
	p10_single.transfer(10, reaction_plate.wells(92).bottom(), reaction_plate.wells(36+i).bottom(0.3), new_tip='never')
	p10_single.blow_out()
p10_single.drop_tip()

# Add 2xMM to 8-part assembly
p10_single.pick_up_tip()
for i in range(12):
	p10_single.transfer(4, reaction_plate.wells(91).bottom(), reaction_plate.wells(48+i).bottom(0.3), new_tip='never')
	p10_single.blow_out()
p10_single.drop_tip()

p10_single.pick_up_tip()
for i in range(12):
	p10_single.transfer(4, reaction_plate.wells(90).bottom(), reaction_plate.wells(60+i).bottom(0.3), new_tip='never')
	p10_single.blow_out()
p10_single.drop_tip()

'''
	# Add 2xMM to 2-part assembly for varies volumes (5,10,20 uL)
	count = 0
	for i in range(3):
	for j in range(8):
	p10_single.pick_up_tip()
	p10_single.transfer(4*(2**i), reaction_plate.wells(89+(3*i)).bottom(0.3), reaction_plate.wells(3*j+count).bottom(0.3), new_tip='never')
	p10_single.blow_out()
	p10_single.drop_tip()
	count += 1
	
	# Add 2xMM to 5-part assembly for varies volumes (5,10,20 uL)
	count = 0
	for i in range(3):
	for j in range(8):
	p10_single.pick_up_tip()
	p10_single.transfer(2.5*(2**i), reaction_plate.wells(88+(3*i)).bottom(0.3), reaction_plate.wells(27+3*j+count).bottom(0.3), new_tip='never')
	p10_single.blow_out()
	p10_single.drop_tip()
	count += 1
	
	# Add 2xMM to 8-part assembly for varies volumes (5,10,20 uL)
	count = 0
	for i in range(3):
	for j in range(8):
	p10_single.pick_up_tip()
	p10_single.transfer(1*(2**i), reaction_plate.wells(87+(3*i)).bottom(0.3), reaction_plate.wells(54+3*j+count).bottom(0.3), new_tip='never')
	p10_single.blow_out()
	p10_single.drop_tip()
	count += 1
	'''

#This function checks the existance of DNA parts and returns for well location of the parts
def find_dna(name, dna_plate_map_dict, dna_plate_dict):
	"""Return a well containing the named DNA."""
	for plate_name, plate_map in dna_plate_map_dict.items():
		for i, row in enumerate(plate_map):
			for j, dna_name in enumerate(row):
				if dna_name == name:
					well_num = 8 * j + i
					return(dna_plate_dict[plate_name].wells(well_num))
	raise ValueError("Could not find dna piece named \"{0}\"".format(name))

#This function checks if the DNA parts exist in the DNA plates and returns for well locaion of output DNA combinations
def find_combination(name, combinations_to_make):
	"""Return a well containing the named combination."""
	for i, combination in enumerate(combinations_to_make):
		if combination["name"] == name:
			return reaction_plate.wells(i)
	raise ValueError("Could not find combination \"{0}\".".format(name))

combinations_by_part = {}
for i in combinations_to_make:
	name = i["name"]
	for j in i["parts"]:
		if j in combinations_by_part.keys():
			combinations_by_part[j].append(name)
		else:
			combinations_by_part[j] = [name]

#This section of the code combines and mix the DNA parts according to the combination list
for part, combinations in combinations_by_part.items():
	part_well = find_dna(part, dna_plate_map_dict, dna_plate_dict)
	combination_wells = [find_combination(x, combinations_to_make) for x in combinations]
	p10_single.pick_up_tip()
	while combination_wells:
		if len(combination_wells) > 5:
			current_wells = combination_wells[0:5]
			combination_wells = combination_wells[5:]
		else:
			current_wells = combination_wells
			combination_wells = []
		p10_single.aspirate(2 * len(current_wells), part_well.bottom(0.5))
		for i in current_wells:
			p10_single.dispense(2, i.bottom(0.5))
		if combination_wells:
			p10_single.mix(2, 10, wash_0.bottom(0.5))
			p10_single.blow_out()
			p10_single.mix(2, 10, wash_1.bottom(0.5))
			p10_single.blow_out()
	# Two washing steps are added to allow recycling of the tips
	p10_single.drop_tip()

# Incubate rxns for 2 hr (moclo), seal the Reaction Plate with adhesive film
'''
	Seal the Reaction Plate to avoid liquid evaporation.
	Discard the empty PCR tubes previously containing buffer, ligase, and restriction enzyme.
	Remove the Input_DNA_Plate from the Deck Space. Remaining DNA may be saved by sealing the Input_DNA_Plate with adhesive film and storing at -20°C
'''

'''start_time = time.time()
temp_deck.set_temperature(37)
p10_single.delay(minutes=120)
num_cols = math.ceil(num_rxns/8.0)

#Adding 4 ul of water halfway through.
	p10_single.pick_up_tip()
	for i in combinations_to_make:
	num_parts = len(i["parts"])
	# Add an extra 4 ul of water for evaporation.
	water_to_add = 4
	well = find_combination(i["name"], combinations_to_make)
	p10_single.transfer(water_to_add, water.bottom(), well.bottom(0.5), new_tip='never')
	p10_single.mix(4, 10, well.bottom(0.5))
	p10_single.mix(2, 10, wash_0.bottom(0.5))
	p10_single.blow_out()
	p10_single.mix(2, 10, wash_1.bottom(0.5))
	p10_single.blow_out()
	p10_single.drop_tip()
	time_elapsed = time.time() - start_time
	p10_single.delay(seconds=(120*60 - time_elapsed))

temp_deck.set_temperature(4)
'''
temp_deck.deactivate()
robot.pause()

# _____________________Start the Transformation protocol________________________________
'''
	The robotic liquid handler would automatically pause when the modular cloning protocol is completed, indicated by the blue light (cooling down) of the Temperature Module.
	Remove the reaction_plate from the top of the Temperature Module and place it on Slot 7 which is assigned to the post_moclo_reaction_plate at the beginning of the protocol, and remove the adhesive film.
	Before beginning the Cell Transformation protocol, load 10 μL/well of chemically competent cells into a new Bio-Rad 96 Well Plate and place it on the top of the Temperature Module.
	Be sure to un-pause the robot after completing all the steps listed above!
	'''

#for i in range(0, num_cols):
#Using letters for rows of custom container to maintain backwards compatibility.
#p300_multi.aspirate(10, reagents_plate.wells('A' + str(i + 1)).bottom(0.5))
#p300_multi.dispense(10, reaction_plate.wells(48 + i*8).bottom(0.5))
#p300_multi.drop_tip()

'''Use 10uL of competent cells with 2uL of each DNA rections from the MoClo protocol'''

# At the beginning of this for loop, the Bio-Rad 96 Well Plate (reaction_plate) now contains 10 uL of competent cells per wells (well 1 through 88)
p10_single.pick_up_tip()
# Add 2 ul of rxns to comp cells
for i in range(0, num_rxns):
	p10_single.transfer(2, post_moclo_reaction_plate.wells(i).bottom(0.5), reaction_plate.wells(i).bottom(0.5), new_tip='never')
	p10_single.mix(4, 10, reaction_plate.wells(i).bottom(0.5))
	p10_single.blow_out()
	p10_single.mix(2, 10, wash_0.bottom(0.5))
	p10_single.blow_out()
	p10_single.mix(2, 10, wash_1.bottom(0.5))
	p10_single.blow_out()
# Two washing steps are added to allow recycling of the tips
p10_single.drop_tip()

# Incubate at 4C, then heat shock.
'''Be sure to un-pause the robot in between the heat shock'''
temp_deck.set_temperature(4)
p10_single.delay(minutes=30)
temp_deck.set_temperature(42)
p10_single.delay(minutes=1)
temp_deck.set_temperature(4)
p10_single.delay(minutes=5)

# Add soc.
p300_multi.pick_up_tip()
for i in range(0, num_cols):
	p300_multi.transfer(150, soc.bottom(), reaction_plate.cols(i).bottom(1), new_tip='never')
	p300_multi.mix(2, 150, reaction_plate.cols(i).bottom(1))
	p300_multi.blow_out()
	p300_multi.mix(2, 300, wash_0.bottom())
	p300_multi.blow_out()
	p300_multi.mix(2, 300, wash_1.bottom())
	p300_multi.blow_out()
# Two washing steps are added to allow recycling of the tips
p300_multi.drop_tip()

# Grow for 1 hr, seal the plate with adhesive film to avoid evaporation
temp_deck.set_temperature(37)
p10_single.delay(minutes=60)
robot.pause()
temp_deck.deactivate()
'''
	Remove the adhesive film from the Reaction Plate before perceeding to Cell Plating.
	Be sure to un-pause the robot after removing the adhesive film!
	'''

# Dilute the recovered transformation reactions and start plating
'''All recovered transformation reactions are diluted to 10% of its original concentration before plating'''
#Dilution
p300_multi.pick_up_tip()
for i in range(0, num_cols):
	p300_multi.transfer(157, reaction_plate.cols(i).bottom(0.5), liquid_waste.bottom(), new_tip='never')
	p300_multi.blow_out()
	p300_multi.mix(2, 300, wash_0.bottom())
	p300_multi.blow_out()
	p300_multi.mix(2, 300, wash_1.bottom())
	p300_multi.blow_out()
	# Two washing steps are added to allow recycling of the tips
	p300_multi.transfer(45, soc.bottom(), reaction_plate.cols(i).bottom(0.5), new_tip='never')
	p300_multi.mix(2, 50, reaction_plate.cols(i).bottom(0.5))
	p300_multi.blow_out()
	p300_multi.mix(2, 300, wash_0.bottom())
	p300_multi.blow_out()
	p300_multi.mix(2, 300, wash_1.bottom())
	p300_multi.blow_out()
# Two washing steps are added to allow recycling of the tips
p300_multi.drop_tip()

#Plating
'''
	Un-comment line 328, and line 333 through line 339 if you're plating for triplicate of each reactions.
	'''
#count2 = 0
for i in range(0, num_cols):
	p300_multi.pick_up_tip()
	p300_multi.transfer(10, reaction_plate.cols(i).bottom(0.5), agar_plate.cols(i).bottom(0.3), new_tip='never')
	#p300_multi.transfer(1, reaction_plate.cols(i).bottom(0.5), agar_plate.cols(i).bottom(-1), new_tip='never')
	#p300_multi.transfer(9, reaction_plate.cols(i).bottom(0.5), agar_plate.cols(i+count2).bottom(0.8), new_tip='never')
	#p300_multi.transfer(1, reaction_plate.cols(i).bottom(0.5), agar_plate.cols(i+count2).bottom(-1), new_tip='never')
	#p300_multi.transfer(9, reaction_plate.cols(i).bottom(0.5), agar_plate.cols(i+count2+1).bottom(0.8), new_tip='never')
	#p300_multi.transfer(1, reaction_plate.cols(i).bottom(0.5), agar_plate.cols(i+count2+1).bottom(-1), new_tip='never')
	#p300_multi.transfer(9, reaction_plate.cols(i).bottom(0.5), agar_plate.cols(i+count2+2).bottom(0.8), new_tip='never')
	#p300_multi.transfer(1, reaction_plate.cols(i).bottom(0.5), agar_plate.cols(i+count2+2).bottom(-1), new_tip='never')
	#count2 += 2
	p300_multi.drop_tip()


'''
	# Dilute and plate.
	def spread_culture(source, dest, soc, dilute_after=True):
	p300_multi.mix(2, 150, source.bottom(0.5))
	p300_multi.aspirate(10, source.bottom(0.5))
	p300_multi.dispense(9, dest.top())
	p300_multi.dispense(1, dest.bottom(-1))
	if dilute_after:
	p300_multi.transfer(120, source.bottom(0.5), liquid_waste.bottom(), new_tip='never')
	p300_multi.mix(2, 300, wash_0.bottom())
	p300_multi.blow_out()
	p300_multi.mix(2, 300, wash_1.bottom())
	p300_multi.blow_out()
	p300_multi.transfer(120, soc, source.bottom(0.5), new_tip='never')
	
	for i in range(0, num_cols):
	agar_plate = agar_plates[i // 3]
	agar_well_num = (i % 3) * 8 * 4
	p300_multi.pick_up_tip()
	source = reaction_plate.wells(48 + i * 8)
	spread_culture(source, agar_plate.wells(agar_well_num), soc)
	spread_culture(source, agar_plate.wells(agar_well_num + 8), soc)
	spread_culture(source, agar_plate.wells(agar_well_num + 16), soc)
	spread_culture(source, agar_plate.wells(agar_well_num + 24), soc, dilute_after=False)
	p300_multi.drop_tip()
	'''





