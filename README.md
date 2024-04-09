# About
This repository has all my gmsh geometries that I made for the Artemis program, to be used for PIC simulations.

# First Time Setup

1. Although using a python virtual environment is not necessary, theres some problems if you install the gmsh library globally (it gets confused when you run gmsh the program) so I recommend creating a virtual environment using
`$ python3 -m venv venv`

2. Activate the virtual environment using
`$ source venv/bin/activate`

3. Install packages (just gmsh in our case)
`$ pip install -r requirements.txt`

4. Run the script with
`$ venv/bin/python ./blue_moon.py`
this assures that we are using the right version of python with the right packages installed

5. To exit the virtual environment when you are done:
`$ deactivate`

# Normal Use

1. Activate the virtual environment
`$ source venv/bin/activate`

2. Run the script
`$ venv/bin/python blue_moon.py`

3. To exit the virtual environment when you are done:
`$ deactivate`