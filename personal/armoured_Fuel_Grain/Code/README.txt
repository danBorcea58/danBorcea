

In order to run the code you shall first install the required libraries.
You can easily do it in one step by pasting inside the project's directory:

python3 -m venv venv
source venv/bin/activate
install pipreqs
pipreqs . --force
pip install -r requirements.txt


gyroidGenerator is the main script. Depending on its input quantities the code might take
from around 10s to more than 1min. The output chart and STLs can be found in Output folder.
The STL file can be then directly imported into any 3D printing slicer.


