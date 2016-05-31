#/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

file="$BASEDIR/../source/revision.cpp"
newrev=$(git rev-parse HEAD)
oldrev=`grep -oP "(?<=\").*?(?=\";)" < $file`
echo "old: $oldrev"
echo "new: $newrev"
if [ "$newrev" != "$oldrev" ] ; then 
  printf "#include \"revision.hpp\"\n\
const std::string kGitRevision = \"$newrev\";\n" > $file
  echo "Revision updated"
else
  echo "Nothing to do"
fi
