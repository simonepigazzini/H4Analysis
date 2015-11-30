#!/bin/bash

queue="1nd"
cfg="cfg/TemplatesMaker.cfg"
exe="bin/TemplatesMaker"

TEMP=`getopt -o q:c:e: --long queue:,cfg:,exe: -n 'submitBatch.sh' -- "$@"`

if [ $? != 0 ] ; then echo "Options are wrong..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

while true; do

  case "$1" in
    -q | --queue ) queue="$2"; shift 2 ;;
    -c | --cfg ) cfg="$2"; shift 2 ;;
    -e | --exe ) exe="$2"; shift 2 ;;
#    -d | --dryrun ) dryrun=1; shift;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done

launchDir=`pwd`

mkdir -p jobs
mkdir -p logs

jobName=`basename ${exe%.*}`"_"`basename ${cfg%.*}`
jobFile="jobs/${jobName}".sh
touch ${jobFile}

cat <<EOF >> ${jobFile}
#!/bin/bash
cd ${launchDir}
source scripts/setup_lxplus.sh
${exe} ${cfg}
EOF

echo "Launching job ${jobName} from dir ${launchDir}"
bsub -q ${queue} -o logs/${jobName}.log -e logs/${jobName}.err -J ${jobName} < ${jobFile}
