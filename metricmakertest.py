"""Unit test for metric maker"""

import metricmaker
import unittest

class MetricTests(unittest.TestCase):
    def testTooManyEmptyTimeCourses(self):
        self.assertTrue(False,'Need resolution on whether time course is going to be allowed to have None in the array')
        time_course = [0,None,None,None,1,2,None,None,None,None ]
        time_points = [1,2,3,4,5,6,7,8,9,10]
        TC_label = None
        sensitivity = None
        metric_select = None
        self.assertRaises(metricmaker.TooManyEmptyValues,
                          metricmaker.generate,
                          time_course, time_points, TC_label, sensitivity, metric_select)
    
    def testValidTimeCourses(self):
        self.assertTrue(False,'Need resolution on whether time course is going to be allowed to have None in the array')
        time_course = [0,1,2,1,1,2,None,None,None,None]
        time_points = [1,2,3,4,5,6,7,8,9,10]
        TC_label = None
        sensitivity = None
        metric_select = None
        metricmaker.generate(time_course, time_points, TC_label, sensitivity, metric_select)

    def testZerosFunction(self):
        z = metricmaker.zeros(1,1)
        self.assertEqual([[0]], z)
        z = metricmaker.zeros(1,2)
        self.assertEqual([[0,0]], z)
        z = metricmaker.zeros(2,2)
        self.assertEqual([[0,0],[0,0]], z)
    
    def testMean(self):
        time_course = [2.219707744, 1.058746906, 1.259369068, 1.400074223, 1.994297498, 1.235476606]
        mean = metricmaker.mean(time_course)
        self.assertEqual(1.5279453408333332, mean)
        
    def testTrapz(self):
        time_points = [1, 2, 3, 4, 5]
        time_course = [10,20,30,40,50]
        value = metricmaker.trapz(time_points, time_course)
        self.assertEqual(120,
                         value)
        #time_course = [2.219707744, 1.058746906, 1.259369068, 1.400074223, 1.994297498]
        time_course = [2.2191, 1.0587, 1.2594, 1.4001, 1.9943]
        self.assertEqual(5.8248999999999995,
                         metricmaker.trapz(time_points, time_course))
        time_points = [0, 0.083333333333333329, 0.25, 0.5, 1, 1.5, 2, 4, 8, 12, 16, 20, 24]
        time_course = [1.2672097275,1.590688972,1.686169412,0.668113898,0.454789704667,0.577257919667,0.747457979667,1.04653912833,1.50108848933,1.58878468067,2.12204150567,1.92499022267,2.41681636767]
        self.assertEqual(38.824680533169918,
                         metricmaker.trapz(time_points, time_course))

    def testDiff(self):
        y = [1.2, 2.5, 0.4, 1, 2]
        self.assertEqual([1.3000, -2.1000, 0.6000, 1.0000],
                         metricmaker.diff(y))
    def testSort(self): 
        data =         [1, 5, 3, 2, 7, 8, 1, 2]
        sorted_data =  [1, 1, 2, 2, 3, 5, 7, 8]
        sorted_index = [0, 6, 3, 7, 2, 1, 4, 5]
        d,i = metricmaker.sort(data)
        # make sure original array wasn't modified
        self.assertEqual([1, 5, 3, 2, 7, 8, 1, 2],
                         data)
        self.assertEqual(sorted_data,
                         d)
        self.assertEqual(sorted_index,
                         i)
        d,i = metricmaker.sort([1, 1, 1, 1, 1, 1])
        self.assertEqual([1, 1, 1, 1, 1, 1],
                         d)
        self.assertEqual([0, 1, 2, 3, 4, 5],
                         i)
                         
    def test_peak_valley(self):
        time_course = [1.2672097275,1.590688972,1.686169412,0.668113898,0.454789704667,0.577257919667,0.747457979667,1.04653912833,1.50108848933,1.58878468067,2.12204150567,1.92499022267,2.41681636767]
        self.assertEqual([-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1],
                         metricmaker.calculate_peak_valley(time_course))
                                          
    def testPeakFinder(self):
        time_points = [0, 0.083333333333333329, 0.25, 0.5, 1, 1.5, 2, 4, 8, 12, 16, 20, 24]
        time_course = [1.2672097275,1.590688972,1.686169412,0.668113898,0.454789704667,0.577257919667,0.747457979667,1.04653912833,1.50108848933,1.58878468067,2.12204150567,1.92499022267,2.41681636767]
        sensitivity = 0.1
        self.assertEqual([-1, 0, 1, 0, -1, 0, 0, 0, 0, 0, 1, -1, 1],
                         metricmaker.peak_finder(time_course, sensitivity))
        # Test at a different sensitivity level
        sensitivity = 0.5
        self.assertEqual([-1, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 1],
                         metricmaker.peak_finder(time_course, sensitivity))
        # Test at a different sensitivity level
        sensitivity = 1.0
        self.assertEqual([-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
                         metricmaker.peak_finder(time_course, sensitivity))
    
    def test_find(self):
        peaks = [1, 0, 1, 0, -1, 0, -1, 0, 1, -1]
        self.assertEqual([0, 2, 8],
                         metricmaker.find(peaks, 1))              
if __name__ == "__main__":
    unittest.main()