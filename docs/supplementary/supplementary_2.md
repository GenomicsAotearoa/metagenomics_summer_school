# NeSI File system, Working directory and Symlinks

The part of the operating system responsible for managing files and directories is called the **file system**. It organizes our data into files, which hold information, and directories (also called ‘folders’), which hold files or other directories.

Directories are like places — at any time while we are using the shell, we are in exactly one place called our current working directory. Commands mostly read and write files in the current working directory, i.e. ‘here’, so knowing where you are before running a command is important.

!!! pied-piper "NeSI Filesystem (For Researchers)"

    All HPC platforms have custom File Systems for the purpose of **general use** and **admin**. 
    
    <center>
    ![image](../theme_images/nesi_filesystem.png){width="450"}
    </center>