#!/bin/bash
#
# This script parses the output of 'quadrule ... +ntries ...' and lists the
# minimum weight of each good rule, as well as its filename and line number.
#
# $Id: quadrule-minwgt.sh,v 1.27 2016/08/26 08:50:07 zlb Exp $

files="$@"
test -z "$files" && files="[stdin]"

for f in $files; do
  if test "$f" = "[stdin]"; then
    f=""
  else
    file -L "$f" 2>/dev/null | grep -qw "text" || continue # skip non text file
  fi
  awk '
    #{printf "=== %d: %s\n", NR, $0 >"/dev/stderr"}	# for debugging

    # This function evaluates a numerical expression in "arg"
    function eval(arg) {
	# return "arg" if it does not contain "<number>[+-*/]<number>"
	if (!match(arg, "[0-9.][-+*/][0-9.]"))
	    return sprintf("%0.16f", arg)
	cmd=sprintf("echo \"%s\" | sed -e \"%s\" | bc -l",
			arg, "s/[EeDd]\\([-+]\\?[0-9]*\\)/*10^(\\1)/g")
	cmd | getline res
	return res
    }

    BEGIN	{flag = 0; w = 1; c = 1; line = 0}
    /^\*    /	{flag = 0; next}		# skip bad rule
    /Dup/ {
	if (!flag) next
	#---------- clean up $0, only leave the number
	sub(".*Dup[1-4]*[(]", ""); sub("[)].*", "")
	D = eval($1)
	if (w > D) w = D
    }
    /Perm|Cons/ {
	if (!flag) next
	pre=$0
	sub(".*Perm.*", "Perm", pre)
	sub(".*Cons.*", "Cons", pre)
	#---------- retrieve "X" from "PermX" or "ConsX", store it in "star"
	star = $0; sub(".*" pre, "", star); sub("[(].*", "", star)
	#---------- clean up $0, only leave the numbers
	sub(".*" pre "[1-4]*[(]", ""); sub("[)].*", ""); gsub(",", " ")
	#---------- process the stars
	if (star == 1111) narg = 3
	else if (star == 111 || star == 211) narg = 2
	else narg = 1
	while (NF < narg) {
	    # add next line
	    if (!getline tmp) {
		printf "Unexpected EOF in \"'$f'\" while processing Perm%s\n",
			star >"/dev/stderr"
		flag = 0;
		exit 1
	    }
	    gsub("[)};]", "", tmp); gsub(",", " ", tmp)
	    $0 = $0 " " tmp
	}
	A = eval($1)
	if (c > A) c = A
	if (star == 4 || star == 3 || star == 2) {
	    # do nothing
	} else if (star == 31) {
	    B = 1.0-3*A
	    B = eval(B)
	    if (c > B) c = B
	} else if (star == 22) {
	    B = 0.5*(1.0-2*A)
	    B = eval(B)
	    if (c > B) c = B
	} else if (star == 211) {
	    B = eval($2)
	    if (c > B) c = B
	    C = 1.0-2*A-B
	    C = eval(C)
	    if (c > C) c = C
	} else if (star == 1111) {
	    B = eval($2)
	    if (c > B) c = B
	    C = eval($3)
	    if (c > C) c = C
	    D = 1.0-A-B-C
	    D = eval(D)
	    if (c > D) c = D
	} else if (star == 21) {
	    B = 1.0-2*A
	    B = eval(B)
	    if (c > B) c = B
	} else if (star == 111) {
	    B = eval($2)
	    if (c > B) c = B
	    C = 1.0-A-B
	    C = eval(C)
	    if (c > C) c = C
	} else if (star == 11) {
	    B = 1.0-A
	    B = eval(B)
	    if (c > B) c = B
	} else {
	    printf "Invalid star (%s:%d): %s\n", "'$f'", NR, star >"/dev/stderr"
	    exit
	}
    }
    /^[23]D [0-9]* point order|^static FLOAT QUAD_[123]D_P[0-9]*_wts/ {
	if (flag)
	    printf "wgt: %0.16f  coor: %0.16f  stars: %s (%s:%d)\n",
			w, c, stars, "'$f'", line
	line = NR; flag = 1; w = 1; c = 1; stars = "n.a."
	if (match($0, "[23]D [0-9]* point order")) {
	    stars = substr($7,2)
	    gsub(",", "", stars)
	}
    }
    END {
	if (flag)
	    printf "wgt: %0.16f  coor: %0.16f  stars: %s (%s:%d)\n",
			w, c, stars, "'$f'", line
    }
  ' $f
done | sort -n -k2,4
