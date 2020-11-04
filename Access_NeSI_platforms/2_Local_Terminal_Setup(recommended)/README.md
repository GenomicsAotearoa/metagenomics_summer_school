

> Following instructions are for Mac,Linux and Windows with WSL enabled

1. In a new local terminal run `mkdir -p ~/.ssh/sockets`

2. Open your *ssh* config file with  `nano ~/.ssh/config` and add the following (copy & paste or type it in). Replace `myusername` with your **NeSI username**

   ```
   Host *
       ControlMaster auto
       ControlPath ~/.ssh/sockets/ssh_mux_%h_%p_%r
       ControlPersist 1
   
   Host ga-vl01
      User myusername
      Hostname ga-vl01.mahuika.nesi.org.nz
      ProxyCommand ssh -W %h:%p lander
      ForwardX11 yes
      ForwardX11Trusted yes
      ServerAliveInterval 300
      ServerAliveCountMax 2
   
   Host lander
      User myusername  
      HostName lander.nesi.org.nz
      ForwardX11 yes
      ForwardX11Trusted yes
      ServerAliveInterval 300
      ServerAliveCountMax 2
   
   Host mahuika
      User myusername  
      Hostname login.mahuika.nesi.org.nz
      ProxyCommand ssh -W %h:%p lander
      ForwardX11 yes
      ForwardX11Trusted yes
      ServerAliveInterval 300
      ServerAliveCountMax 2
   ```

   

3. Close and save by pressing `Ctrl+X`  >  `Y` > `Enter`

4. Ensure the permissions are correct by running `chmod 600 ~/.ssh/config`