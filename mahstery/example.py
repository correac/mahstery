from eagle_util import readEAGLE
from mahstery import run

print('Example: Calculating best-fit expression from EAGLE DMO L0100N1504 simulation')

#! Please, indicate the path of EAGLE database for mahstery
# e.g. '/Users/Camila/Downloads/'
dir = '/Users/Camila/Dropbox/mahstery/'

print('Reading EAGLE data, wait a bit..')
Mz, z, c = readEAGLE(dir)

print('Running mahstery')
run(Mz, z, c)
print('Done. mahstery output mahstery_output.png file')
