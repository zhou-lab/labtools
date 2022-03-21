# sf means "Sync From"

# Ex1: remote > local, because /mnt/isilon is a server path
# from /mnt/isilon/zhoulab/labprojects/a/b/c to ~/zhoulab/labprojects/a/b/c
# sf /mnt/isilon/zhoulab/labprojects/a/b/c

# Ex2: remote > local, because /home/zhouw3 is a server path
# from /mnt/isilon/zhoulab/labprojects/a/b/c to ~/zhoulab/labprojects/a/b/c
# sf /home/zhouw3/zhoulab/labprojects/a/b/c

# Ex3: remote > local, because of explicit r:
# from ~/zhoulab/labprojects/a/b/c to ~/zhoulab/labprojects/a/b/c
# sf r:~/zhoulab/labprojects/a/b/c

# Ex4: local > remote, because /Users/zhouw3 is a local path
# from /User/zhouw3/zhoulab/labprojects/a/b/c to /home/zhouw3/zhoulab/labprojects/a/b/c
# sf /Users/zhouw3/zhoulab/labprojects/a/b/c

# Ex5: local > remote, default if there is no way to tell direction
# from ~/zhoulab/labprojects/a/b/c to ~/zhoulab/labprojects/a/b/c
# sf ~/zhoulab/labprojects/a/b/c

### INSTALLATION
### insert the following to your .zshrc/.bashrc
### ----------------
# LOCAL_HOME="/Users/zhouw3"
# REMOTE_HOME="/home/zhouw3"
# HPC_NAME="hpc5"
# source <(curl -s https://raw.githubusercontent.com/zhou-lab/labwiki/master/Snippets/20210326_sync_HPC_data.sh)
### ----------------

function sf() {

  ## setting default if not given
  [[ -z "$LOCAL_HOME" ]] && LOCAL_HOME="/Users/zhouw3"
  [[ -z "$REMOTE_HOME" ]] && REMOTE_HOME="/home/zhouw3"
  [[ -z "$REMOTE_HOME2" ]] && REMOTE_HOME2="/mnt/isilon"
  [[ -z "$REMOTE_SCR" ]] && REMOTE_SCR="/scr1/users/zhouw3/"
  [[ -z "$HPC_NAME" ]] && HPC_NAME="hpc"

  from=$1
  if [[ $from =~ ^$LOCAL_HOME ]]; then # from local to remote
    to=${from/$LOCAL_HOME/$HPC_NAME":"$REMOTE_HOME}
  elif [[ $from =~ ^$REMOTE_HOME2 ]]; then # from remote to local
    to=${from/$REMOTE_HOME2/$LOCAL_HOME}
    from=$HPC_NAME":"$from
  elif [[ $from =~ ^$REMOTE_SCR ]]; then 
    to=${from/"$REMOTE_SCR"/"$LOCAL_HOME/scr1_zhouw3/"}
    from=$HPC_NAME":"$from
  elif [[ $from =~ ^$REMOTE_HOME ]]; then # from remote to local
    to=${from/$REMOTE_HOME/$LOCAL_HOME}
    from=$HPC_NAME":"$from
  elif [[ $from =~ ^r: ]]; then # from remote to local
    from=${from/r:/$HPC_NAME":"}
    from=${from/\~/$REMOTE_HOME}
    from=${from/$LOCAL_HOME/$REMOTE_HOME}
    to=${from/$HPC_NAME":"/}
    to=${to/#\~/$LOCAL_HOME}
    to=${to/$REMOTE_HOME/$LOCAL_HOME}
  elif [[ $from =~ ^~ ]]; then # from local to remote
    from=${from/\~/$LOCAL_HOME}
    to=${to/#\~/$REMOTE_HOME}
    to=$HPC_NAME":"$from
  else
    return 1;
  fi

  echo "From: "$from
  echo "To:   "$to

  if [[ $from =~ ^$LOCAL_HOME ]]; then
    ssh $HPC_NAME mkdir -p $(dirname ${to/$hpc_name":"/})
    rsync -av --update $from $to
  elif [[ $from =~ ^$HPC_NAME ]]; then
    mkdir -p $(dirname $to)
    rsync -av --update $from $to
  fi
}
