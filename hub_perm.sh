#!/bin/bash
H=/expanse/projects/sebat1/jsebat/cnv_bmi_autistic_traits
cd "$H" || exit 1
chgrp -R jsebat-group "$H"; chmod -R g+rX,o-rwx "$H"; find "$H" -type d -exec chmod g+s {} +
for s in inputs derived code results; do echo "$s files: $(find "$H/$s" -type f | wc -l)  size: $(du -sh "$H/$s" | cut -f1)"; done
ls -la "$H"
