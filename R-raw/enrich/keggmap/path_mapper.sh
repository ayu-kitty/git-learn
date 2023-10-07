#!/bin/bash
# by The Coder, 20220926

#参数1：差异文件名
#参数2：anno文件名
#参数3：keggmap路径
#参数4：颜色

xlsmap(){
	path=$(basename "$1" |sed 's/.xls$//'|sed 's/\r//')
	html=$km/$path.html
	png=$km/$path.png
	HTML="$dirn"/$path.html
	all=`echo $col | cut -d "," -f 3`
	up=`echo $col | cut -d "," -f 1`
	down=`echo $col | cut -d "," -f 2`
	
	if [ -s $html ] && [ -s $png ]; 
	then
	  cp $png "$dirn"/
	else
		echo "WARNING: missing $path!"
		rm "$1"
		continue
	fi
	awk 'NR==FNR{if($0~"^<head>") p1=FNR; if($0~"<!-- pathway image start -->") p2=NR; if($0~"<!-- pathway image end -->") p3=NR;if($0~"<!-- pathway image url end -->") p4=NR;next}
    (FNR<=p1 ||( FNR>=p2&&FNR<=p3)||FNR>p4){print $0}' $html $html |sed 's#"/dbget-bin/#"http://www.genome.jp/dbget-bin/#g' |
	sed 's#"/kegg-bin/#"http://www.kegg.jp/kegg-bin/#' |sed "s:./${path}_files/$path.png:$path.png:"|
    sed "s:/kegg/pathway/.*.png\" srcset=:$path.png\" srcset=:"|sed "s:/kegg/pathway/.*.png:$path.png:;s/ title=/\ttitle=/" |sed 's:src="./:src=":'|sed 's/1x"/"/'|
	awk 'BEGIN{FS=OFS="title="} NR==FNR{a[$2]=1; next} 
	(($1~"circle" ||$1~"rect" || $1~"poly") && $2!~":"){for(i in a){ p=match($2, "[^[:alnum:]]"i"[^[:alnum:]]"); if(p>0) $0=$0"title="i}}{print $0}' "$1" - | 
	awk 'BEGIN{FS=OFS="title="}NR==FNR{a[$2]=(a[$2] ? a[$2]""$1 : $2""$1); next}(NF>2){$2="title=\""a[$3]; for(i=4; i<=NF; i++){ $2=$2"&#13;"a[$i] };$2=$2"\" />"; $0=$1$2} {print $0}' "$1" - | 
	sed "s:./$path.png 1x:$path.png:" > "$HTML"

	awk 'BEGIN{FS=OFS="coords=|data-coords|title="} ($1~"circle" ||$1~"rect" || $1~"poly"){print $2"\t"match($4, "Up")"\t"match($4, "Down")}' "$HTML" | 
	sed 's/"//g'| 
	awk 'BEGIN{FS=OFS="\t"} ($2+$3)>0{if($2*$3>0) $4="'$all'"; if($3==0) $4="'$up'";if($2==0) $4="'$down'"; print $1,$4}' | 
	awk 'BEGIN{FS=OFS="\t"} {split($1, a, ",")} 
	length(a)==4{print "region", $2 , "46x18+"a[1]"+"a[2] }
	length(a)==3{print "circ", $2 ,a[1]" "a[2]","a[1]+4" "a[2] }
	length(a)==8{print "draw", $2 , "line;"a[1]","a[2]";"a[3]","a[4] }' |
	while read a b c; do
    if [ "$a" == "region" ]; then
        convert "$dirn"/$path.png -quiet -auto-level -fuzz 80% -fill $b -region $c -opaque white "$dirn"/$path.png
    fi
    if [ "$a" == "draw" ]; then
        draw=$(echo $c | sed 's/;/ /g')
        convert "$dirn"/$path.png -quiet -auto-level -stroke $b -strokewidth 2 -draw "$draw"  "$dirn"/$path.png
    fi
    if [ "$a" == "circ" ]; then
        convert "$dirn"/$path.png -quiet -auto-level -fill $b -draw 'circle '"$c"'' "$dirn"/$path.png
    fi
	done

	rm "$1"

}
export -f xlsmap
km=$3
export km
dirn="$dirn"
export dirn
col=$4
export col
keggmap1(){
	dirn=$(basename "$diff" | sed 's/-diff-.*//')
	mkdir "$dirn"

	awk 'BEGIN{FS=OFS="\t"} {p=match($2, "."); if(p>0) $2=substr($2,1,p+4)} 
	NR>1{print $1, "FC: "$2, $3}' "$diff" | 
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2"\t"$3; next} 
	(a[$1]){print $0,a[$1]}' $anno - | 
	sed 's/path://g; s/PATH://g' | 
	awk -v dirn="$dirn" 'BEGIN{FS="\t"; OFS="    "} {split($5, a, ",")} 
	{for(i in a) print "&#13;",$1,$2,$3"title="$4 > dirn"/"a[i]".xls"}'
	
  if [ -e "$dirn"-list.xls ]
  then
    ls $dirn | grep -v -w -f "$dirn"-list.xls | xargs -I {} rm -f $dirn/'{}'
    rm "$dirn"-list.xls
  fi
	ls "$dirn"/*.xls |/opt/conda/bin/parallel -j 10 xlsmap

}
export -f keggmap1

diff="$1"
export diff
anno=$2
export anno

keggmap1 "$1" "$2" "$3" "$4"
rm "$1" "$2"