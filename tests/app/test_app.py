import unittest
import os
from tests.config.definitions import ROOT_DIR
from app.app import App
from sdk.moveapps_io import MoveAppsIo
import pandas as pd
import random

print("Conducting unit tests. ")

class MyTestCase(unittest.TestCase):

    def setUp(self) -> None:
        os.environ['APP_ARTIFACTS_DIR'] = os.path.join('C:/Users/matth/hi_outdoor_recreation_analysis', 'tests/resources/output')
        self.sut = App(moveapps_io=MoveAppsIo())  

    #Checks the input file compared to the output after it was run
    def test_app_returns_input(self):
        # prepare
        expected = pd.read_pickle("./resources/samples/martesJul24.pickle")

        #config: dict = {}
        config: dict = {"keys": "Key-Pair-Id=APKAIDPUN4QMG7VUQPSA&Policy=eyJTdGF0ZW1lbnQiOiBbeyJSZXNvdXJjZSI6Imh0dHBzOi8vaGVhdG1hcC1leHRlcm5hbC0qLnN0cmF2YS5jb20vKiIsIkNvbmRpdGlvbiI6eyJEYXRlTGVzc1RoYW4iOnsiQVdTOkVwb2NoVGltZSI6MTcyMjc5ODQyOX0sIkRhdGVHcmVhdGVyVGhhbiI6eyJBV1M6RXBvY2hUaW1lIjoxNzIxNTc0NDI5fX19XX0_&Signature=iC1sIRxZz3Azaalz~0wPCCzxmXOSEH3hBpfMDOqucqZl7QyAA9UreV74YBQCJJ~fDeLAAqFBpcHb0UTdo50rxXfiWw5S0XB~Hbtbk9DOm5OzpjtqlgDa9ro8tmWxZx~QJSvhXB1YuawWyBc0seeOsZMBlpn30IVCeLaN-QwrRRnGs15irxLubqGwuk7WMdyIgVnYnMQ3qdB-eWpSZ1iez11taiI5TL~TBdgsZ3pCVRa4RKO4kVqj859o406NnKZLgnslz--lI9204wquMk9i8NXDMtfBfs2Qk8eXj-afUAUhBIFrlUnxkkAmDhAHP2kJZZX83fPMBW-scb42QMPwCg__"}

        # execute
        actual = self.sut.execute(data=expected, config=config)

        #verify
        expectedList = []
        actualList = []
        expectedLength = len(expected.to_point_gdf())
        actualLength = len(actual.to_point_gdf())
        expectedGdf = expected.to_point_gdf()
        actualGdf = actual.to_point_gdf()
        actCols = len(actualGdf.iloc[0])
        expectedList.append(expectedGdf.iloc[expectedLength-1].iloc[0])
        actualList.append(actualGdf.iloc[actualLength-1].iloc[0])

        #checks a 500 column row coordinates at random within the dataset excluding the 7 columns that I have changed or added inlcuding point, distance, intensity, geometry, band 1, band 2, band 3, and band 4
        for i in range(500):
            randRow = random.randint(0, actualLength - 1)
            randCol = random.randint(0, actCols - 8)
            if pd.isna(actualGdf.iloc[randRow].iloc[randCol]):
                continue
            else:
                expectedList.append(expectedGdf.iloc[randRow].iloc[randCol])
                actualList.append(actualGdf.iloc[randRow].iloc[randCol]) 

        #check the list
        self.assertListEqual(expectedList, actualList)

        #geometries should be unequal between the expected and actual
        self.assertNotEqual(expectedGdf.iloc[actualLength-1].iloc[actCols-7], actualGdf.iloc[actualLength-1].iloc[actCols-7])
        
if __name__ == '__main__':
    unittest.main()