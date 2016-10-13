Installation instructions for Optizelle

1.  Drag the folder ${CPACK_PACKAGE_LONG_NAME} to Applications

2.  Set the environment variables for desktop use of Optizelle.  Copy the file
    called

        /Applications/${CPACK_PACKAGE_LONG_NAME}/share/optizelle/com.optimojoe.optizelle.plist

    to either 

        ~/Library/LaunchAgents/

    which enables Optizelle for a local user or

        /Library/LaunchAgents/

    which enables Optizelle for all users.  Then, reboot the computer.  Note,
    by default, the folder ~/Library is hidden inside of Finder.  To access it,
    use the menu option Go->Go to Folder..., type ~/Library, and then click the
    Go button.

3.  (Optional) Set the environment variables for use of Optizelle over an ssh
    session.  For bash, add

        export PYTHONPATH=$PYTHONPATH:/Applications/${CPACK_PACKAGE_LONG_NAME}/share/optizelle/python
        export MATLABPATH=$MATLABPATH:/Applications/${CPACK_PACKAGE_LONG_NAME}/share/optizelle/matlab
        export OCTAVE_PATH=$OCTAVE_PATH:/Applications/${CPACK_PACKAGE_LONG_NAME}/share/optizelle/octave

    to ~/.bash_login
