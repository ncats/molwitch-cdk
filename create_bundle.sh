#!/bin/sh 
echo "Please enter the name for the bundle file:" 
read filename 
echo "You entered: $filename" 
mkdir bundle 
cp pom.xml bundle 
cp target/*.jar bundle 
cd bundle 
directory=$(pwd) 
for file in "$directory"/*; do 
    if [ -f "$file" ]; then 
        # Perform an action on each file 
        echo "Signing file: $file" 
        gpg -ab $file         
    fi 
done 
jar cvf $filename.jar *
echo "$filename.jar is created." 