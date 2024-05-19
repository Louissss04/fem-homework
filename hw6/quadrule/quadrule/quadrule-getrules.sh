#!/bin/bash
#
# This script lists all good rules in the output of 'quadrule ... +ntries ...'.
#
# $Id: quadrule-getrules.sh,v 1.6 2016/07/25 07:00:07 zlb Exp $

files="$@"
test -z "$files" && files="[stdin]"

for f in $files; do
  if test "$f" = "[stdin]"; then
    f=""
  else
    file -L "$f" 2>/dev/null | grep -qw "text" || continue # skip non text files
  fi
  awk '
    BEGIN	{flag = 0; line = 0}
    /^[23]D .* point order/ {
	if (flag) print rule
	line = NR; flag = 1
	rule=sprintf("%s (%s:%s)", $0, "'$f'", line)
	next
    }
    END 	{if (flag) print rule}
    /^\*    /	{flag = 0; next}		# bad rule
    /^ *Dup/ 	{if (flag) rule=rule "\n" $0; next}
    /^ *Perm[0-9]/ {if (flag) rule=rule "\n" $0; next}
    /^ *Cons[0-9]/ {if (flag) rule=rule "\n" $0; next}
  ' $f
done
