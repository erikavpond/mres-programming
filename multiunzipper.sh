#should unzip multiple files into folder matching original name, and return relevant error message.

for zipfile in *.zip
do
  dirname=`echo $zipfile | sed 's/\.zip$//'`
  if mkdir "$dirname"
  then
    if cd "$dirname"
    then
      unzip ../"$zip"
      cd ..
    else
      echo "Could not unpack $zip - cd failed"
    fi
  else
    echo "Could not unpack $zip - mkdir failed"
  fi
done
