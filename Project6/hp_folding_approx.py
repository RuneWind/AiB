
def calculate_matches(hp_str):

	even_s = [i for i in range(len(hp_str)) if i%2 == 0 and hp_str[i] == "h"]
	odd_s = [i for i in range(len(hp_str)) if i%2 == 1 and hp_str[i] == "h"]
	#print(even_s, len(even_s))
	#print(odd_s, len(odd_s))
	
	even_odd_score = 0
	odd_even_score = 0
	
	odd_s_rev = odd_s[::-1]
	even_s_rev = even_s[::-1]

	# Match evens with reverse odds
	i = 0
	even = even_s[i]
	odd = odd_s_rev[i]

	while even < odd:
		even_odd_score += 1
		
		i += 1
		even = even_s[i]
		odd = odd_s_rev[i]

	# Match reverse evens with odds
	i = 0
	even = even_s_rev[i]
	odd = odd_s[i]

	while odd < even:
		odd_even_score += 1
		
		i += 1
		even = even_s_rev[i]
		odd = odd_s[i]

	if even_odd_score > odd_even_score:
		return zip(even_s[0: even_odd_score], odd_s_rev[0: even_odd_score])

	if even_odd_score < odd_even_score:
		return list(zip(odd_s[0: odd_even_score], even_s_rev[0: odd_even_score]))


#def calculate_free_energy():









hp_str = "hhhhhhhhhhhhphphpphhpphhpphpphhpphhpphpphhpphhpphphphhhhhhhhhhhh"
print(calculate_matches(hp_str))

