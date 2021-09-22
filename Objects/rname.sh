for f in *.pickle; do cp "$f" "$(echo "$f" | sed s/.pickle/_lowpt.pickle/)"; done
