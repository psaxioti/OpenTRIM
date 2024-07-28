#! /usr/bin/env sh

# Usage
#   mingw64-deploy.sh install_folder distname
#
# input:
#   install_folder : the folder with the iradina++ install files.
#                    This is typically [build_folder]/iradinapp
#   distname : e.g., iradinapp-0.1.6

# the script does the following
#   1. Copy exe and dll to a folder named "distname"
#   1. find all /ucrt64 dlls that exe and dll depend on
#   2. copy these to distname folder
#   3. create a zip archive named distname 

INSTALLPATH=$1
DISTNAME=$2
CURFLDR=$PWD

mkdir $DISTNAME
cp $INSTALLPATH/bin/*.* $DISTNAME
cp $INSTALLPATH/lib/*.dll $DISTNAME

cd $DISTNAME 

printf "ldd iradina++.exe\n"
list=$(ldd ./iradina++.exe | sed 's/[^\/]*\(\/[^ ]*\)/\1\n/' | grep ucrt64)
for dll in $list;
do
  dll_lst="$dll_lst $dll"
done

printf "ldd ./libiradinapp.dll\n"
list=$(ldd ./libiradinapp.dll | sed 's/[^\/]*\(\/[^ ]*\)/\1\n/' | grep ucrt64)
for dll in $list;
do
  dll_lst="$dll_lst $dll"
done

printf "ldd ./libiondedx.dll\n"
list=$(ldd ./libiondedx.dll | sed 's/[^\/]*\(\/[^ ]*\)/\1\n/' | grep ucrt64)
for dll in $list;
do
  dll_lst="$dll_lst $dll"
done

printf "ldd ./libxs_zbl.dll\n"
list=$(ldd ./libxs_zbl.dll | sed 's/[^\/]*\(\/[^ ]*\)/\1\n/' | grep ucrt64)
for dll in $list;
do
  dll_lst="$dll_lst $dll"
done

printf "ldd ./libxs_lj.dll\n"
list=$(ldd ./libxs_lj.dll | sed 's/[^\/]*\(\/[^ ]*\)/\1\n/' | grep ucrt64)
for dll in $list;
do
  dll_lst="$dll_lst $dll"
done

printf "ldd ./libxs_krc.dll\n"
list=$(ldd ./libxs_krc.dll | sed 's/[^\/]*\(\/[^ ]*\)/\1\n/' | grep ucrt64)
for dll in $list;
do
  dll_lst="$dll_lst $dll"
done

printf "ldd ./libxs_moliere.dll\n"
list=$(ldd ./libxs_moliere.dll | sed 's/[^\/]*\(\/[^ ]*\)/\1\n/' | grep ucrt64)
for dll in $list;
do
  dll_lst="$dll_lst $dll"
done

# remove duplicates
dll_lst=`echo $dll_lst | tr ' ' '\n' | sort | uniq`

for dll in $dll_lst;
do
  cp $dll .
done

cd $CURFLDR

tar -a -cf $DISTNAME.zip $DISTNAME

