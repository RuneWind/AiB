import os
import time

running_times = []

hp_strs = ["hhppppphhppphppphp", "hphphhhppphhhhpphh", "phpphphhhphhphhhhh", "hphpphhphpphphhpphph", "hhhpphphphpphphphpph", "hhpphpphpphpphpphpphpphh", "pphpphhpppphhpppphhpppphh", "ppphhpphhppppphhhhhhhpphhpppphhpphpp", "pphpphhpphhppppphhhhhhhhhhpppppphhpphhpphpphhhhh", "hhphphphphhhhphppphppphpppphppphppphphhhhphphphphh", "pphhhphhhhhhhhppphhhhhhhhhhphppphhhhhhhhhhhhpppphhhhhhphhphp", "hhhhhhhhhhhhphphpphhpphhpphpphhpphhpphpphhpphhpphphphhhhhhhhhhhh", "hhhhpppphhhhhhhhhhhhpppppphhhhhhhhhhhhppphhhhhhhhhhhhppphhhhhhhhhhhhppphpphhpphhpphph", "pppppphphhppppphhhphhhhhphhpppphhpphhphhhhhphhhhhhhhhhphhphhhhhhhppppppppppphhhhhhhpphphhhpppppphphh", "ppphhpphhhhpphhhphhphhphhhhpppppppphhhhhhpphhhhhhppppppppphphhphhhhhhhhhhhpphhhphhphpphphhhpppppphhh"]
for i in range(len(hp_strs)):
	print(i+1, "\t", end="")
	
	t_start = time.time()
	match_indecies = calculate_matches(hp_strs[i])
	fold_string = create_fold_string(match_indecies, hp_strs[i])
	t_end = time.time()
	running_times.append(t_end - t_start)