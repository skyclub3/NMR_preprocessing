# AWK library
# Add the option -i to include this library when you use awk

# check rotamer states
# input: number from -180 to 180; output: gp (gauche +), gm (gauche -), tr (trans)
function state(str){rstr=(str>0&&str<120?"gp":(str<0&&str>-120?"gm":"tr"));return rstr;}

# trim string ; remove the characters in front or back
function trim(str){sub(/^[ \t]+/,"",str);sub(/[ \t]+$/,"",str);return str}
function trimparen(str){sub(/^[(]+/,"",str);sub(/[)]+$/,"",str);return str}
function trimquot(str){sub(/^['"]+/,"",str);sub(/['"]+$/,"",str);return str}
function trimslash(str){sub(/^[/]+/,"",str);sub(/[/]+$/,"",str);return str}
#it must consider that it removes the quots in the middle ; the next version removes quots in the first and end 
function remquot(str){gsub(/"/,"",str);return str}
#062320 add trim white spaces, \n \t \r
function trimwhite(str){sub(/^[ \t\r\n]+/,"",str);sub(/[ \t\r\n]+$/,"",str);return str}
#121520 add new trim removing special characters '")}] ; not working
#function trimspecial(str){sub(/['")}\]]+/,"",str);return str}

#011421 trim brakets, {, }, [, and ], in first and last
function trimbraket(str){sub(/^[\[\{]+/,"",str);sub(/[\]\}]+$/,"",str);return str}

# file exist
function fexist(file){if((getline<file)> 0){close(file);return 1}else{return 0}}
function dexist(dir){if(system("ls -d "dir">&/dev/null")==0){close("ls -d "dir">&/dev/null");return 1}else{return 0}}

# pdbformat (#col=18)
# atom#: 2 ; atomtype: 4 ; alternate: 5 ; residue: 6 ; chain: 8 ; resid: 9 ; insertion: 10 ; x: 12 ; y: 13 ; z: 14 
function pdbformat(){FIELDWIDTHS=" 6 5 1 4 1 3 1 1 4 1 3 8 8 8 6 6 6 4";printpdbformat="%-6s%5s%1s%4s%1s%3s%1s%1s%4s%1s%3s%8s%8s%8s%6s%6s%6s%4s\n"}

# vector arithmetic
function innerprod(a,b){return a[1]*b[1]+a[2]*b[2]+a[3]*b[3]}
function vabs(a){return sqrt(a[1]*a[1]+a[2]*a[2]+a[3]*a[3])}
function vdist(a,b){return sqrt((a[1]-b[1])**2.+(a[2]-b[2])**2.+(a[3]-b[3])**2.)}
# component input
function innerprodcom(ax,ay,az,bx,by,bz){return ax*bx+ay*by+az*bz}
function vabscom(ax,ay,az){return sqrt(ax*ax+ay*ay+az*az)}
function vdistcom(ax,ay,az,bx,by,bz){return sqrt((ax-bx)**2.+(ay-by)**2.+(az-bz)**2.)}
function vectorcom(ax,ay,az,bx,by,bz,retvec){retvec[1]=(bx-ax);retvec[2]=(by-ay);retvec[3]=(bz-az)}

#trigonometric functions
function asin(x) { return atan2(x, sqrt(1-x*x)) }
function acos(x) { return atan2(sqrt(1-x*x), x) }
function atan(x) { return atan2(x,1) }


# if a number has beyond the assigned number of digit, the function removes decimals: sigintnum(a,b) ; a: number, b: a assigned number of digit
# ex) sigintnum(a,3) : number a has digit above 3 (>=3), it has the decimal to be removed.
function sigintnum(a,b){if(length(int(a))>=b){return int(a)}else{return a}}


#functions for amino acid handling
#to use the dictionary r1tor3[], r3tor1[], first run aminoacid(str) ; if str is "help", return the usage  
function aminoacid(str){
  if(str=="help"){
    printf "r1tor3[A]=ALA r3tor1[ALA]=A\n" # help message : it is not working. why? -> worked
    exit
  }else{
    resinfo="/home/jinhyuk/bin/resinfo" #residue info (col3: one letter code ; col2: three letter code)
    while((getline < resinfo)>0){
      r1tor3[$3]=$2; # change one letter code of amino acid to three letter code (it is global variable although it defines in subroutine, aminoacid())
      r3tor1[$2]=$3; # inversely change
    }
  }
}


#math functions
function abs(x){if(x<0){return -x}else{return x}}


#052120 max/min functions ; input array, such as a[]
#default awk (version 4.0.2) cannot include lib.awk by -i option ; use /opt/utils/gawk-5.1.0/gawk which is aliased gawk and awk
# variable names, max , index are used in awk? don't use user variables
function maxarray(ar,    j,result,findj){
  #returns: max pos 
  result=-99999
  for(j in ar){
    if(ar[j]+0==ar[j]&&ar[j]>result){result=ar[j];findj=j}else{} # ar[j]+0==ar[j] is used for filter out non-numbers
  }
  return sprintf("%s %s",result,findj)
}
function minarray(ar,    j,result,findj){
  result=99999
  for(j in ar){
    if(ar[j]+0==ar[j]&&ar[j]<result){result=ar[j];findj=j}else{}
  }
  return sprintf("%s %s",result,findj)
}

#string maninuplation functions
#word(str,n) words(str) ; these are same in gnuplot internal functions ; delimiter: default(blank, space)
function word(str,n,    a){
  na=split(str,a," ")
  return a[n]
}
function words(str,    a){
  na=split(str,a," ")
  return na
}
#general word() words() => gword(str,n,sep) gwords(str,sep)
function gword(str,n,seperator,    a){
  na=split(str,a,seperator)
  return a[n]
}
function gwords(str,seperator){
  na=split(str,a,seperator,    a)
  return na
}

#052520
function str2array(str,array,seperator){
  return split(str,array,seperator)
}

function array2str(array,seperator,    i,j,str){
  str="";
  j=1;
  for(i in array){
    #str=sprintf("%s%c%s",str,(j==1?"":seperator),array[i])
    #102820 unseen chracter (^) is existing ; it can be seen in an editor, Emacs ; because of the character, the variable is not saved (set a = ``) ; so it changed to return the blank, " " or change %c to %s 
    str=sprintf("%s%s%s",str,(j==1?"":seperator),array[i])

    j++
  }
  return str
}

#return true(1) or false(0) of exact matching string in array
function arraymatch(string,array,    x,y){
  for(x in array){
    y[array[x]]
  }
  return string in y
}

#052920, 060120 ; source from https://github.com/e36freak/awk-libs/blob/master/math.awk ; the comments are excerpt from 052520 Google Calendar and the name of 'awk'
#function isnumber(a){return (a~/^[0-9]+$/)} ; 1 2
#function isnumber(a){return (a~/^[-+]?[0-9]+$/)} ; +1 -2   
#function isnumber(a){return (a~/^[-+]?[0-9.]+$/)} ;  +1 -2 +1.0 -1.2
#function isnumber(a){return (a~/^[-+]?[0-9.]+$/ && a !~ /\..*\./)} ; +1 -2 +1.0 -1.2 ;; 8.2.3 (is not number ; this routine is correct but how???)
#function   isnumber(a){return (a+0==a)} ; all ok but only two -inf and -infinity pretends to be number not string ; why????
# see: https://unix.stackexchange.com/questions/584421/inf-in-awk-not-working-the-way-inf-does/584442#584442  
#   awk 'BEGIN{print inf}' => nothing happen ; awk 'BEGIN{print -inf}' ==> 0 (defined) ????
#   or awk 'BEGIN{print +inf}' ==> 0 (defined)
#   awk 'BEGIN{a=-inf;b=inf;print typeof(a),typeof(b)}' ==> number unassigned(untyped)
#    see https://unix.stackexchange.com/questions/584421/inf-in-awk-not-working-the-way-inf-does/584442#584442  
#awk handle +INF, +Infinity, -INF, -Infinity number 0?
#                     not INF, Infinity (bug!!!!!)
#
#gnuplot also handle all of +INF, +Infinity, INF, Infinity, -INF, -Infinity ; number ; with them it fails to stat ; if they are any strings not the upper, it skip the row and do stats

function isnumber(a){
  return (a~/^[-+]?[0-9.]+$/ && a !~ /\..*\./)}

#060420 RMSE, R(Pearson correlation)
function RMSE(arraya,arrayb,    k,l,i,sum){
  k=0 # if it is not initialized here, it enter here by parameter ; for it to work as a local variables, it must initialize 
  l=0
  sum=0
  for(i in arraya) k++
  for(i in arrayb) l++
  if(k!=l){
      printf("Array sizes are not matched")
   return 0
  }else{
  }
  for(i=1;i<=k;i++){
    sum=sum+(arraya[i]-arrayb[i])**2
  }
  return sqrt(sum/k)
}

#060520 transpose of array2d with nrow and ncol ; size: array2d[nrow,ncol]
#return: 2d newarray2d
# ==> by using function is not appropriate ; we cannot accept and printout whole ; it's better to use awk script (by using -f)
#function transpose(array2d,nrow,ncol,newarray2d,    i,j){
#  for(i=1;i<=ncol;i++){
#    for(j=1;j<=nrow;j++){
#      newarray2d[i,j]=sprintf("%s%c",array2d[j,i],(j==nrow?"\n":" "))
#    }
#  }
#}


#062520 Pearson correlation (R)
#PearsonCorr(array1,array2,nrow,    i,sx,sy,ax,ay,xy,xx,yy) ; array1 and 2 is array with number index from 1 to nrow
#see 062520 google calendar ; array index must numbers
#usage: awk -i ~/bin/lib.awk '{array1[NR]=$1;array2[NR]=$2}END{nrow=NR;print PearsonCorr(array1,array2,nrow)}' test.txt ; test.txt has two columns for which we measure R 
function PearsonCorr(array1,array2,nrow,   i,sx,sy,ax,ay,xy,xx,yy){
    sx=0
    sy=0
    for(i=1;i<=nrow;i++){
	sx+=array1[i]
	sy+=array2[i]
    }
    ax=sx/nrow
    ay=sy/nrow
    xy=0
    xx=0
    yy=0
    for(i=1;i<=nrow;i++){
	xy+=(array1[i]-ax)*(array2[i]-ay)
	xx+=(array1[i]-ax)**2
	yy+=(array2[i]-ay)**2
    }
    return xy/sqrt(xx)/sqrt(yy)
}

#version 2
#PearsonCorr2(array1,array2,   i,ii,sx,sy,ax,ay,xy,xx,yy)  ; array1 and 2 are same index (string or number) ; it requires the size of the array
# array index could strings
#usage: awk -i ~/bin/lib.awk '{array1[NR"X"]=$1;array2[NR"X"]=$2}END{nrow=NR;print PearsonCorr2(array1,array2)}' test.txt ; test.txt has two columns for which we measure R
function PearsonCorr2(array1,array2,   i,ii,sx,sy,ax,ay,xy,xx,yy){
    sx=0
    sy=0
    ii=0
    for(i in array1){
	sx+=array1[i]
	sy+=array2[i]
	ii++
    }
    ax=sx/ii
    ay=sy/ii
    xy=0
    xx=0
    yy=0
    for(i in array1){
	xy+=(array1[i]-ax)*(array2[i]-ay)
	xx+=(array1[i]-ax)**2
	yy+=(array2[i]-ay)**2
    }
    return xy/sqrt(xx)/sqrt(yy)
}

#070920 HTML URL encoding function ; change the special characters to HTML URL encoding
# ref: https://www.w3schools.com/tags/ref_urlencode.ASP
function HTMLencoding(str,    changestr,na,a,i){
    na=split(str,a,"")
    changestr=""
    for(i=1;i<=na;i++){
	changestr=sprintf("%s%s",changestr,chr2UTF[a[i]])
    }
    return changestr
}


#define constants? and filename
BEGIN{
    #constants defined
    pi=4.*atan2(1.0,1.0);

    #070920 ; dictionary for character to UTF8
    char2UTF8="/data2/yycho/Protein/yycho/opt/src_packages/char2UTF8.txt" # character to UTF8 (HTML encoding)
    if("z" in chr2UTF){
    }else{
	while((getline < char2UTF8)>0){
	    chr2UTF[substr($0,1,1)]=substr($0,2)
	}
    }

}

#091120 sort string in ascending order and return in array 
function sortstr(string,array,order,   n,a,i,j){
    n=split(string,a," ")
    #see the gawk pdf p320 
    if(order=="asc"){
	PROCINFO["sorted_in"]="@val_num_asc"
    }else if(order=="desc"){
	PROCINFO["sorted_in"]="@val_num_desc"
    }else{
        #usage: awk -i lib.awk 'BEGIN{str="2 3 5";sortstr(str,array,"asc");for(i in array){print i,array[i]}}'
	printf("Usage: sortstr(str,array,order): str: string which has numbers ; array: return sorted array ; order: asc or desc\n")
	exit
    }
    i=1
    for(j in a){
	array[i]=a[j]
	i++
    }
}

#091420 sort array in order and return in array
function sortarray(array1,array2,index2,order,   i,j ){
    if(order=="asc"){
	PROCINFO["sorted_in"]="@val_num_asc"
    }else if(order=="desc"){
	PROCINFO["sorted_in"]="@val_num_desc"
    }else{
	#Usage: awk -i ~/bin/lib.awk 'BEGIN{a[3]=1;a[2]=2;a[1]=3;sortarray(a,array,ind,"asc");for(i in array){print i,ind[i],array[i]}}'
	printf("Usage: sortarray(array1,array2,order): array1: input array which has numbers with index ; array2: return sorted array and index array (index2) ; order: asc or desc\n")
	exit
    }
    i=1
    for(j in array1){
	index2[i]=j
	array2[i]=array1[j]
	i++
    }
}

#012121 add rotamer accuracy function
function rotaccu(chif,chis,tol,     c1,c2,d,s,f){
    c1=(chif<=0.0?chif+360:chif);
    c2=(chis<=0.0?chis+360:chis);
    d=(c1-c2>0?c1-c2:c2-c1);
    d=(d>180?360-d:d);
    return (d<=tol?1:0)
}

#082621 array length
function alen(a,     i,k) {
    k=0
    for (i in a) k++;
    return k;
}

