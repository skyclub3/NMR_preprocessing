BEGIN{
#        Width = 10; 
#        Interval = 1.0;
        Height = 20; 
        Epsilon = 0.000001; 


        MinData =  999999999999; 
        MaxData = -999999999999; 
        MinCount =  999999999999; 
        MaxCount = -999999999999; 
        NumIndex = 0; 
} 
{ 
        # 
        # Pass in only one field of data 
        # The data must be integer, if continuos it will be converted 
        # to integer. 
        # 
        # The data must range greater than 1 because it is converted 
        # to an integer. For data < 1 multiply by a factor to make it 
        # > than 1 
        # 
        Value = $1; 

        # 
        # Get the Max/Min value of the data 
        # 
        if(Value < MinData) { 
                MinData = Value; 
        } 
        if(Value > MaxData) { 
                MaxData = Value; 
        } 

        # 
        # Convert the data value into an int 
        # 
#        ValueInt = int(Value); 
        ValueInt = Value; 


        # 
        # Create a Index of different data values 
        # 
        if(Count[ValueInt] == "") { 
                Index[NumIndex] = ValueInt; 
                NumIndex++; 
        } 


        # 
        # Count the Value in the Histogram Count array 
        # 
        Count[ValueInt]++; 


} 
END { 
	if(MaxData>=0.0){
	  MaxData=int(MaxData)+1;
	}
	if(MaxData<0.0){
	  MaxData=int(MaxData);
	}
	if(MinData>=0.0){
	  MinData=int(MinData);
	}
	if(MinData<0.0){
	  MinData=int(MinData)-1;
	}

	Width=(MaxData-MinData)/Interval;

#        print MaxData,MinData,Width;

        # 
        # Sort the Index array based on its value 
        # 
        QuickSort(Index, 0, NumIndex-1); 

        # 
        # Bin up the data into the reduced scale histogram plot 
        # That is re-scale it 
        # 
        for(i=0;i<Width;i++) { 
                HistX[i] = 0; 
        } 
        Scale = Width/(MaxData-MinData+Epsilon); 
        for(i=0;i<NumIndex;i++) { 
                Col=int((Index[i] - MinData) * Scale); 
                HistX[Col] = HistX[Col] + Count[Index[i]]; 
        } 

        # 
        # Get the Max/Min of the Histogram Count 
        # 
        for(i=0;i<Width;i++) { 
                if(HistX[i] < MinCount) { 
                        MinCount = HistX[i]; 
                } 
                if(HistX[i] > MaxCount) { 
                        MaxCount = HistX[i]; 
                } 
        } 

#Print out
        for(i=0;i<Width;i++) {
#	  Col = MinData + i * (MaxData-MinData)/ (Width-1);
	  Col = MinData + i * Interval
	  print Col,HistX[i];
	}

#        print MinData, MaxData, Width

} 
# 
# Function QuickSort 
# 
function QuickSort(Index, left, right) { 
        # 
        # Do nothing if array contains 
        # fewer than two elements 
        # 
        if(left >= right) { 
                return; 
        } 

        # 
        # Move Partition element to Index[0] 
        # 
        Middle = int((left + right)/2); 
        Swap(Index, left, Middle); 
        last = left; 


        # 
        # Partition 
        # 
        for(i=left+1; i<=right; i++) { 
                if(Index[i] < Index[left]) { 
                        Swap(Index, ++last, i); 
                } 
        } 


        # 
        # Restore partition element 
        # 
        Swap(Index, left, last); 
        QuickSort(Index, left, last-1); 
        QuickSort(Index, last+1, right); 
} 
# 
# Function Swap 
# 
function Swap(Index, i, j) { 
        temp = Index[i]; 
        Index[i] = Index[j]; 
        Index[j] = temp; 
} 
