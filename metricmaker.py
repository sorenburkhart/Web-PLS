"""Convert to and from Roman numerals

 Functions which creates a new metrics based on the characteristics of a timecourse passed to the function
 Translated to python from code written in MATLAB by Jason Kelly, August 2003

 Desciption:
 metric_maker accepts a timecourse and returns a number of metrics based on this timecourse, including
 mean, AUC, max values, decay rates, equilibrium, and activiation slopes.

"""

__author__ = "Soren Burkhart (soren.burkhart@gmail.com)"
__version__ = "$Revision: 0.2 $"
__date__ = "$Date: 2010/02/01 19:40:22 $"
__copyright__ = "Copyright (c) 2010 Soren Burkhart"
__license__ = "Python"

import math
import copy

#Define exceptions
class MetricMakerError(Exception): pass
class TooManyEmptyValues(MetricMakerError): pass
class TimePointsMismatchWithTimeCourse(MetricMakerError): pass
class UnknownError(MetricMakerError): pass

def find(array, value):
    "Searches through array for value. Returns an array of the indexes found"
    return map(None,(i for i in xrange(len(array)) if array[i] == value))
    
def sort(data):
    "Sorts an array.  Returns a copy of the sorted array and the index of the sort order"
    zipped_data = zip(data, xrange(len(data)))
    zipped_data.sort()
    sorted_data = map(None, (v[0] for v in zipped_data))
    sorted_index = map(None,(v[1] for v in zipped_data))
    
    return sorted_data, sorted_index
    
def diff(data):
    "Calculates Difference and approximate derivative."
    diff = []
    for i in range(1, len(data)):
        diff.append(data[i] - data[i-1])
    return diff
    
def trapz(time_points, time_course):
    "Calculates the area under the curve using the trapazoid rule"
    area = 0.0
    for i in range(1,len(time_points)):
        height = time_points[i]-time_points[i-1]
        base = (time_course[i]+time_course[i-1])/2
        area += height * base
        #print "i = %f height = %f base = %f area = %f" % (i,height,base,area)
    return area
    
def mean(values):
    "Returns the mean of a list."
    return sum(values)/len(values)
    
def zeros(rows, columns):
    "Provides a zero matrix with the specified rows and columns"
    matrix = []
    for row in xrange(rows):
        row = []
        for col in xrange(columns):
            row.append(0)
        matrix.append(row)
    return matrix
    
def generate(time_course, time_points, TC_label, sensitivity, metric_select):
    """
     INPUT: One timecourse.
     Input format: metric_maker(time_course, time_points, TC_label, sensitivity, metric_select)

     time_course:
     An array of length N where N = number of timepoints in the timecourse.
     The values in the array correspond to the levels of the assayed variable (i.e. kinase
     activity or portein level).  Format of the array should be:
     Pr_t1 Pr_t2 Pr_t3 ...
     where Pr_t1 is protein at time point #1, #2, etc.

     time_points:
     An array of length N where N = number of timepoints in the timecourse.
     This is simply a listing of the timepoints which correspond to the values
     given in time_course.

     sensitivity:
     A number between 0-1 which determines how large a peak needs to be in order to be considered significant. 
     Peak height is determined by comparing the tip of the peak to the lowest valleys on either side.  Both
     distances must exceed the distance X=sensitivity * (max point in timecourse - min point in timecourse).

     TC_label:
     A string which provides the name of the timecourse measured  (i.e. AKT, JNK)

     metric_select:
     Array of length J where J=number of the metrics supported by metric maker.  (Currently J=8)
     A metric is given a user-defined value of 1 if it is to be included in the metric generation
     and a value of 0 if it is to be omitted.
     current metric_select order: 
     [mean, AUC_whole, Max, Equilibrium, Derivative, Data_Points, AUC_peak, ActivationSLope_peak, DecayRate_peak]
    """
    
    # first check if there are more empty values than measurments
    # if there are too many empty values, then metric_maker exits
    if sum(1 for v in time_course if v == None) >= len(time_points)/2:
        raise TooManyEmptyValues, "found too many empty values in the time course."
    
    # number of timepoints in the timecourse.
    num_timepoints = len(time_course)
    
    #----Initialize New Metric Variables----------------------
    avg = 0;        # averages across the timecourse
    peak = 0;       # peak value of variable across the time course
    AUC = 0;        # itegral of the timecourse 
    maxpeaktime = 0;# time until peaks
    equil = 0;      # equilibrium activity based on the last 25% of timepoints
    deriv = zeros(1,num_timepoints-1); # derivative at each timepoint
    peaktime = 0;   # time of the max point in the timecourse.
    
    # Call the peakfinder to build the peak_valley_course
    peak_valley_course = peak_finder(time_course,sensitivity)

    # Create indexes of the locations of peaks and valleys in the timecourse.
    peak_index = find(peak_valley_course == 1);
    num_of_peaks = length(peak_index);
    valley_index = find(peak_valley_course == -1);
    num_of_valleys = length(valley_index);

def calculate_peak_valley(y_data):
    "Determines the peaks for a timeseries"
    # number of timepoints in the data series
    num_timepoints = len(y_data);
    
    # Initialization of the output array. 
    peak_valley = zeros(1,num_timepoints)[0];
    
    # find all peaks (both insignificant and significant)

    # Differentiate the data matrix at each point in the y_data array, 
    diff_y_data = diff(y_data);

    # set the last_sign var to be equal to the slope of the first point.
    if diff_y_data[0]<0:
        last_sign = -1
    else:
        last_sign = 1
    
    #print "diff_y_data = %s" % diff_y_data
    
    # find the peaks in the timecourse wherever a slope changes from positive to negative (ignoring the endpoints for now)
    for i in range(1,num_timepoints-2):
        if diff_y_data[i]<0:
            sign = -1
        else:
            sign = 1
        
        if (sign == -last_sign):
            #if going from positive to negative slope
            if(sign < last_sign):
                peak_valley[i] = 1
            last_sign = sign
        
    # determine if endpoints are peaks or valleys, the current algorithm considers endpoints peaks or 
    # valleys regardless of whether they are significant.
    if diff_y_data[0] > 0: 
        #if positive slope between first 2 points, 1st point is a valley
        peak_valley[0] = -1
    else:
        #else its a peak
        peak_valley[0] = 1
    
    if diff_y_data[num_timepoints-2]>0: 
        #if positive slope between last 2 points, last point is a peak
        peak_valley[num_timepoints-1] = 1
    else:
        #else its a valley
        peak_valley[num_timepoints-1] = -1
    
    return peak_valley
    
def peak_finder(y_data, sensitivity):
    """
       Function which determines the "significant" peaks/valleys for an input timecourse.
       Translated to python from code written in MATLAB by Jason Kelly, August 2003.

       Description 
       This function excepts a timeseries and a sensitivity level and generates an array 
       which represents the peaks and valleys of the timeseries.

       INPUT: a timeseries and a sensitivity level
       Input format: (y_data, sensitivity)

       y_data:
       An array of values representing a timeseries.  (The actual timepoints of the timeseries
       are not required as they do not come into play in determining the signifiance of peaks in 
       current algorithm.

       sensitivity:
       Value bewteen 0-1 which represents the level of significance required for a peak to be 
       considered significant.  
       Sensitivity * (max value in time series - min value in timeseries) = minimum size for a significant peak.
       "size" of the peak is determined by the distance between the peak and the nearest 2 valleys.

       OUTPUT: array of length M where M equals the number of elements in the y_data array.
       Each position in the array holds a valut of -1, 0 , or 1 and corresponds to a timepoint
       in the y_data array.
       -1 = valley
        0 = insiginificant 
        1 = peak
    """
    # number of timepoints in the data series
    num_timepoints = len(y_data);
    #print "num_timepoints = %d" % num_timepoints
    # the size of the minimum "jump" from a valley to a peak based on the signiificance.
    min_jump = (max(y_data)-min(y_data))*sensitivity;
    #print "min_jump = %f" % min_jump
    
    # Create a combined matrix containing the data points along with their respecive peak/valley characteristic
    # Sorted according to time points
    peak_valley = calculate_peak_valley(y_data)
    combo_sort_by_TP = [y_data, peak_valley]
    #print "Combo sort = %s" % combo_sort_by_TP
        
    sorted_timecourse, sort_TC_index = sort(combo_sort_by_TP[0]);
        
    #print "sorted_timecourse = %s" % sorted_timecourse
    
    # Sort combined matrix by amplitude (time course)
    combo_sort_by_TC = zip(y_data, peak_valley)
    combo_sort_by_TC.sort()
    combo_sort_by_TC.reverse()
    
    sort_TC_index.reverse()
    
    #print        "combo_sort_by_TC = %s\nsort_TC_index = %s" % (combo_sort_by_TC,sort_TC_index)
    
    # Go through all time points
    for i in range(0,num_timepoints):
        #print "NEXT TIMEPOINT i = %d" %i
        #print "combo_sort_by_TP\n%s\ncombo_sort_TC\n%s\nsort_TC_index\n%s" % (combo_sort_by_TP, combo_sort_by_TC, sort_TC_index)
        
        # If the current time point is a peak.
        if combo_sort_by_TC[i][1] == 1:
            j=1
            left_sig=0 #1 = significant ; -1 = insignificant
            # while left side not significant and not the beginning or end points and whose left point isn't at time position 0 (out of bounds)
            
            while (left_sig!=1)&(sort_TC_index[i]!=num_timepoints-1)&(sort_TC_index[i]!=0)&((sort_TC_index[i]-j)!=-1):
                #----------COMPARE TO THE LEFT-------------------------------
                position_in_time = (sort_TC_index[i]-j)
                current_peak = [combo_sort_by_TP[0][sort_TC_index[i]],combo_sort_by_TP[1][sort_TC_index[i]]]
                compare_left = [combo_sort_by_TP[0][sort_TC_index[i]-j],combo_sort_by_TP[1][sort_TC_index[i]-j]]
                #print "compare left"
                #print "position_in_time=%d\ncurrent_peak=%s\ncompare_left=%s\n" % (position_in_time, current_peak, compare_left)
                # if the current peak is greater than the significant jump size from the compare point, it's significant.
                if ((current_peak[0] - compare_left[0])>min_jump):
                    left_sig = 1
                
                # if the comparing point is a peak and is lower then the current peak, then the comparing point is not a significant peak.
                # these rules dont apply to the very first point.
                if ((compare_left[1] == 1)&((sort_TC_index[i]-j)!=1)):
                    if ((current_peak[0]-compare_left[0]))>0:
                        combo_sort_by_TP[1][sort_TC_index[i]-j]=0;
                    else: #if the comparing peak is higher than the current peak, then the current peak is not a significant peak.
                        left_sig=-1;
                        combo_sort_by_TP[1][sort_TC_index[i]]=0;

                # if it reaches the beginning of the time course and the beginning is a peak, then current peak is insignificant.
                # if beginning is a valley than current peak is significant.
                if ((sort_TC_index[i]-j)==0):
                    if (compare_left[1] == 1):
                        left_sig=-1
                        combo_sort_by_TP[1][sort_TC_index[i]]=0
                     
                    if (compare_left[1] == -1):
                        left_sig = 1
                     
                    if (compare_left[1] == 0):
                         raise UnknownError, 'problem: shouldnt get here --  end peaks were not set properly. (Must be either a peak or valley)'  

                j=j+1
                #print 'while loop left i=%d num_timepoints=%d left_sig=%d sort_TC_index[i]=%d sort_TC_index[i]-j=%d' % (i,num_timepoints,left_sig,sort_TC_index[i],sort_TC_index[i]-j)

            #print 'end compare left while'
            j=1
            right_sig=0 #1 = significant ; -1 = insignificant

            #print 'before right i=%d num_timepoints=%d right_sig=%d sort_TC_index[i]=%d sort_TC_index[i]-j=%d' % (i,num_timepoints,right_sig,sort_TC_index[i],sort_TC_index[i]-j)
            #--------------COMPARE TO THE RIGHT-------------------------------
            while (right_sig!=1)&(sort_TC_index[i]!=num_timepoints-1)&(sort_TC_index[i]!=0)&((sort_TC_index[i]+j)<=num_timepoints-1):
                position_in_time = (sort_TC_index[i]+j)
                #print 'num_timepoints = %d\ni = %d\nj = %d\n\nsort_TC_index[i] = %d\ncombo_sort_by_TP%s' % (num_timepoints,i,j,sort_TC_index[i],combo_sort_by_TP)
                current_peak = [combo_sort_by_TP[0][sort_TC_index[i]],combo_sort_by_TP[1][sort_TC_index[i]]]
                compare_right = [combo_sort_by_TP[0][sort_TC_index[i]+j],combo_sort_by_TP[1][sort_TC_index[i]+j]]
                
                #print "compare right"
                #print "position_in_time=%d\ncurrent_peak=%s\ncompare_right=%s\n" % (position_in_time, current_peak, compare_right)
                
                # if the current peak is greater than the significant jump size from the compare point, its sugnificant.
                if ((current_peak[0] - compare_right[0])>min_jump):
                    right_sig = 1
                
                # if the comparing point is a peak and is lower then the current peak, then the comparing point is not a significant peak
                # these rules dont apply to the very last point.
                #print "compare_right[1] = %d sort_TC_index[i]+j != num_timepoints -1  %d != %d" % (compare_right[1],sort_TC_index[i]+j,num_timepoints-1)
                if ((compare_right[1] == 1)&((sort_TC_index[i]+j)!=num_timepoints-1)):
                    #print "INSIDE CHECKPOINT"
                    if ((current_peak[0]-compare_right[0])>0):
                        #print "peak"
                        combo_sort_by_TP[1][sort_TC_index[i]+j]=0
                    else:
                        #print "NO PEAK"
                        # if the comparing peak is higher than the current peak, then the current peak is not a significant peak
                        right_sig = -1
                        combo_sort_by_TP[1][sort_TC_index[i]]=0
                #print "SKIP CHECKPOINT"
                # if it reaches the end of the time course and the end is a peak, then current peak is insig
                # if end is a valley than current peak is significant
                if ((sort_TC_index[i]+j)==num_timepoints-1):
                    if (compare_right[1] == 1):
                        right_sig=-1;
                        combo_sort_by_TP[1][sort_TC_index[i]]=0
                    
                    if (compare_right[1] == -1):
                        right_sig = 1

                    if (compare_right[1] == 0):
                        raise UnknownError, 'problem: shouldnt get here --  end peaks were not set properly. (Must be either a peak or valley)'  
                
                j=j+1

            # Rebuild the combo sorted by time course to include the lost peaks (if peak was significant it wont be lost)
            #[sorted_timecourse,sort_TC_index] = sort(combo_sort_by_TP(:,1));
            #% sort by amplitude (time course)
            #combo_sort_by_TC = combo_sort_by_TP(sort_TC_index,:);
            #combo_sort_by_TC = combo_sort_by_TC((end:-1:1),:);
            #sort_TC_index = sort_TC_index((end:-1:1),:);
            sorted_timecourse, sort_TC_index = sort(combo_sort_by_TP[0]);
        
            #print "sorted_timecourse = %s" % sorted_timecourse
    
            # Sort combined matrix by amplitude (time course)
            combo_sort_by_TC = zip(y_data, peak_valley)
            combo_sort_by_TC.sort()
            combo_sort_by_TC.reverse()
    
            sort_TC_index.reverse()
    #print "combo_sort_by_TP %s" % combo_sort_by_TP
    #---------SET VALLEYS--------------------------------
    # Valleys are set as the lowest point between two peaks.

    # there hasn't been a first peak yet
    no_first = 1 
    # initialize the low valley as the highest point, so valley will be guaranteed to be lower.
    low_valley = max(combo_sort_by_TP[0])
    #print 'low_valley = %f' % low_valley
    for i in range(0,num_timepoints):
        if (combo_sort_by_TP[1][i] == 1): 
            if (no_first == 1):
                no_first = 0
                low_valley = max(combo_sort_by_TP[0])
            else:
                # if its a peak (but not the first one) set the low valley for the previous section and reset the low valley
                combo_sort_by_TP[1][low_valley_position] = -1
                low_valley = max(combo_sort_by_TP[0])
        if (combo_sort_by_TP[0][i]<low_valley):
            low_valley = combo_sort_by_TP[0][i]
            low_valley_position = i
    
    #print        "combo_sort_by_TP %s" % combo_sort_by_TP
    return combo_sort_by_TP[1]
    