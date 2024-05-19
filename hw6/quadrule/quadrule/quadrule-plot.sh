#!/bin/bash
#

echo 1>&2 "Reading rule from stdin:"
awk '/),$/{print; next} /,$/{printf "%s", $0; next} {print}' | \
	grep -i 'Perm\|Cons' | sed -e 's/%\|\*//g' -e 's/),/)/g' | \
	sed -e 's/[Pp][Ee][Rr][Mm]\|[Cc][Oo][Nn][Ss]/\\Perm/g' | \
	sed -e 's/[	 ]//g' | \
	tee /tmp/quadrule-plot-tmp.tex
cat <<END >/tmp/quadrule-plot-tmp-main.tex
\\def\\inputfile{/tmp/quadrule-plot-tmp.tex}\\input quadrule-plot.tex
END
pdflatex /tmp/quadrule-plot-tmp-main.tex
mv -v quadrule-plot-tmp-main.pdf quadrule-plot.pdf
/bin/rm -f /tmp/quadrule-plot-tmp.tex /tmp/quadrule-plot-tmp-main.tex \
	quadrule-plot-tmp-main.*
