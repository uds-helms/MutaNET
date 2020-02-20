graph [
	directed 1
	graphics [
	]
	LabelGraphics [
	]
	node [
		id 3
		label "gn6 (3,1,1)"
		category_of_interest "MDR efflux regulator"
	]
	node [
		id 4
		label "lt7 (1,1,1)"
		category_of_interest "-"
	]
	node [
		id 5
		label "lt8 (2,1,1)"
		category_of_interest "antibiotic"
	]
	node [
		id 13
		label "gn16 (2,1,0)"
		category_of_interest "MDR efflux pump"
	]
	node [
		id 14
		label "gn16_1 (2,1,0)"
		category_of_interest "MDR efflux pump"
	]
	node [
		id 17
		label "gn20 (0,1,0)"
		category_of_interest "-"
	]
	node [
		id 20
		label "gn26 (0,0,0)"
		category_of_interest "antibiotic"
	]
	node [
		id 21
		label "gn27 (2,2,0)"
		category_of_interest "-"
	]
	edge [
		id 1
		source 21
		target 14
		label "?"
	]
	edge [
		id 2
		source 3
		target 13
		label "A"
	]
	edge [
		id 3
		source 21
		target 13
		label "?"
	]
	edge [
		id 4
		source 3
		target 14
		label "A"
	]
	edge [
		id 5
		source 3
		target 5
		label "R"
	]
	edge [
		id 6
		source 5
		target 4
		label "O"
	]
	edge [
		id 7
		source 4
		target 3
		label "O"
	]
	edge [
		id 8
		source 17
		target 5
		label "?"
	]
	edge [
		id 9
		source 17
		target 3
		label "?"
	]
]