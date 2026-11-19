#!/bin/bash
H=/expanse/projects/sebat1/jsebat/cnv_bmi_autistic_traits; cd "$H"
ok=0; bad=0
while IFS=$'\t' read -r set file bytes md5; do md5=${md5%$''};  [ "$set" = "set" ] && continue
  f="$H/$set/$file"; if [ -f "$f" ] && [ "$(md5sum "$f" | cut -d' ' -f1)" = "$md5" ]; then ok=$((ok+1)); else bad=$((bad+1)); echo "MISMATCH $set/$file"; fi
done < DATA_MANIFEST.tsv
echo "verified $ok, mismatched $bad"
