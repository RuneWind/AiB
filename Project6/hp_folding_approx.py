import os
import time

def calculate_matches(hp_str):

	# Lists of even indices and odd indices of h's in input "hp"-string
	even_s = [i for i in range(len(hp_str)) if i%2 == 0 and hp_str[i] == "h"]
	odd_s = [i for i in range(len(hp_str)) if i%2 == 1 and hp_str[i] == "h"]
	
	# Lists in reverse lists
	odd_s_rev = odd_s[::-1]
	even_s_rev = even_s[::-1]
	
	# Number of matchings when taking the evens from left and the odds from right
	even_odd_score = 0
	# Number of matchings when taking the odds from left and the evens from right
	odd_even_score = 0
	
	# Match evens from left with odds from right, i.e. match evens with reverse odds
	i = 0
	while even_s[i] < odd_s_rev[i]:
		even_odd_score += 1
		i += 1
		if i >= min(len(odd_s_rev), len(even_s)):
			break

	# Match odds from left with evens from right, i.e. match odds with reverse evens
	i = 0
	while odd_s[i] < even_s_rev[i]:
		odd_even_score += 1
		i += 1
		if i >= min(len(odd_s), len(even_s_rev)):
			break
		
	# Return tuples of matchings in list from the largest number of matchings
	if even_odd_score > odd_even_score:
		return list(zip(even_s[0: even_odd_score], odd_s_rev[0: even_odd_score]))

	if even_odd_score <= odd_even_score:
		return list(zip(odd_s[0: odd_even_score], even_s_rev[0: odd_even_score]))


def create_fold_string(match_indecies, hp_str):
	free_energy = 0
	fold_str = ""

	# Add e's corresponding to bases before first match
	fold_str += match_indecies[0][0] * "e"

	min_overhang = min(match_indecies[0][0], len(hp_str)-match_indecies[0][1]-1)
	for i in range(1, min_overhang+1):
		if hp_str[match_indecies[0][0]-i] == "h" and hp_str[match_indecies[0][1]+i] == "h":
			free_energy += 1


	# Make left to right portion of fold string
	for i in range(0, len(match_indecies)-1):
		# In case of loop
		if match_indecies[i+1][0] - match_indecies[i][0] > 2:
			fold_str += "n" * ((match_indecies[i+1][0] - match_indecies[i][0] - 1) // 2)
			fold_str += "e"
			fold_str += "s" * ((match_indecies[i+1][0] - match_indecies[i][0] - 1) // 2)
			fold_str += "e"
			if hp_str[match_indecies[i+1][0]-1] == "h":
				free_energy += 1
				if hp_str[match_indecies[i][1]-1] == "h":
					free_energy += 1
			free_energy += 1


		# In case of base pairs
		else:
			fold_str += 2 * "e"
			if hp_str[match_indecies[i][0]+1] == "h" and hp_str[match_indecies[i][1]-1] == "h":
				free_energy += 1
			free_energy += 1

	tail_length = ((match_indecies[len(match_indecies)-1][1] - match_indecies[len(match_indecies)-1][0] - 1) // 2)

	for i in range(tail_length):
		upper = hp_str[match_indecies[len(match_indecies)-1][0] + i]
		lower = hp_str[match_indecies[len(match_indecies)-1][1] - i]
		if upper == "h" and lower == "h":
			free_energy += 1

	# Add tail of structure to fold string
	fold_str += tail_length * "e"
	fold_str += "s"
	fold_str += tail_length * "w"


	# Make right to left portion of fold string
	for i in range(len(match_indecies)-1, 0, -1):
		# in case of loop
		if match_indecies[i-1][1] - match_indecies[i][1] > 2:
			fold_str += "s" * ((match_indecies[i-1][1] - match_indecies[i][1] - 1) // 2)
			fold_str += "w"
			fold_str += "n" * ((match_indecies[i-1][1] - match_indecies[i][1] - 1) // 2)
			fold_str += "w"
			if hp_str[match_indecies[i-1][1]-1] == "h":
				free_energy += 1

		# In case of base pairs
		else:
			fold_str += 2 * "w"

	fold_str += (len(hp_str) - match_indecies[0][1] - 1) * "w"
	
	return fold_str, -free_energy




# Measure running time of implementation
for i in range(5):
	running_times = []
	lens = []
	
	hp_strs = ["hhppppphhppphppphp", "hphphhhppphhhhpphh", "phpphphhhphhphhhhh", "hphpphhphpphphhpphph", "hhhpphphphpphphphpph", "hhpphpphpphpphpphpphpphh", "pphpphhpppphhpppphhpppphh", "ppphhpphhppppphhhhhhhpphhpppphhpphpp", "pphpphhpphhppppphhhhhhhhhhpppppphhpphhpphpphhhhh", "hhphphphphhhhphppphppphpppphppphppphphhhhphphphphh", "pphhhphhhhhhhhppphhhhhhhhhhphppphhhhhhhhhhhhpppphhhhhhphhphp", "hhhhhhhhhhhhphphpphhpphhpphpphhpphhpphpphhpphhpphphphhhhhhhhhhhh", "hhhhpppphhhhhhhhhhhhpppppphhhhhhhhhhhhppphhhhhhhhhhhhppphhhhhhhhhhhhppphpphhpphhpphph", "pppppphphhppppphhhphhhhhphhpppphhpphhphhhhhphhhhhhhhhhphhphhhhhhhppppppppppphhhhhhhpphphhhpppppphphh", "ppphhpphhhhpphhhphhphhphhhhpppppppphhhhhhpphhhhhhppppppppphphhphhhhhhhhhhhpphhhphhphpphphhhpppppphhh"]
	for i in range(len(hp_strs)):
		t_start = time.time()
		match_indecies = calculate_matches(hp_strs[i])
		fold_string = create_fold_string(match_indecies, hp_strs[i])[0]
		t_end = time.time()
		running_times.append(t_end - t_start)
		lens.append(len(hp_strs[i]))
	
	print(running_times)
print(lens)
	

#match_indecies  = calculate_matches(hp_str)
##print(match_indecies)
#fold_string = create_fold_string(match_indecies, hp_str)
#
#print(len(hp_str), len(fold_string))
#print(match_indecies)
#print(fold_string)
#
#os.system("python hpview3k.py " + hp_str + " " + fold_string )


