for i in $(seq -f "%3g" 0 200); do
#echo "checking lxplus ${i} "
ssh -o ConnectTimeout=1 \
 -o PreferredAuthentications=gssapi-with-mic \
 -o GSSAPIAuthentication=yes \
 -o StrictHostKeyChecking=no \
 -o LogLevel=quiet \
 lxplus9$i.cern.ch 'output=$(tmux list-sessions);
hostname
if [[ "${output}" == *"no server running on"* ]]; then
    echo $output
else
    echo "Killing the server"
    tmp=$(tmux kill-server)
    
fi
'
done
