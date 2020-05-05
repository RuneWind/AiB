
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

	# Match odds from left with evens from right, i.e. match odds with reverse evens
    i = 0
    while odd_s[i] < even_s_rev[i]:
        odd_even_score += 1
        i += 1
        
    # Return tuples of matchings in list from the largest number of matchings
    if even_odd_score > odd_even_score:
        return zip(even_s[0: even_odd_score], odd_s_rev[0: even_odd_score])

    if even_odd_score <= odd_even_score:
        return list(zip(odd_s[0: odd_even_score], even_s_rev[0: odd_even_score]))


#def calculate_free_energy(match_indecies):
	




hp_str = "hhhhhhhhhhhhphphpphhpphhpphpphhpphhpphpphhpphhpphphphhhhhhhhhhhh"
match_indecies  = calculate_matches(hp_str)
print(match_indecies)
#calculate_free_energy(match_indecies)




