
.PHONY: all
all : convexHull435

convexHull435 : convexHull435.o
	g++ convexHull435.o -o convexHull435

convexHull435.o : convexHull435.cpp
	g++ -std=c++17 -c convexHull435.cpp

ConvexHull.o : ConvexHull.cpp
	g++ -std=c++17 -c ConvexHull.cpp

.PHONY: clean
clean :
	@rm -f *.o convexHull435 hull_G* hull_J* hull_Q*
	@rm -f *10*.txt

#generates all hull_J_*.txt files with Jarvis March.
.PHONY: runJ
runJ : convexHull435
	time -o output/Jarvis/J_test_rectangle10.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J rectangle10.txt
	time -o output/Jarvis/J_test_rectangle100.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J rectangle100.txt
	time -o output/Jarvis/J_test_rectangle1000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J rectangle1000.txt
	time -o output/Jarvis/J_test_rectangle10000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J rectangle10000.txt
	time -o output/Jarvis/J_test_rectangle100000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J rectangle100000.txt
	time -o output/Jarvis/J_test_rectangle1000000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J rectangle1000000.txt
	time -o output/Jarvis/J_test_circle10.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J circle10.txt
	time -o output/Jarvis/J_test_circle100.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J circle100.txt
	time -o output/Jarvis/J_test_circle1000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J circle1000.txt
	time -o output/Jarvis/J_test_circle10000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J circle10000.txt
	time -o output/Jarvis/J_test_circle100000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J circle100000.txt
	time -o output/Jarvis/J_test_circle1000000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J circle1000000.txt
	time -o output/Jarvis/J_test_onCircle10.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J onCircle10.txt
	time -o output/Jarvis/J_test_onCircle100.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J onCircle100.txt
	time -o output/Jarvis/J_test_onCircle1000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J onCircle1000.txt
	time -o output/Jarvis/J_test_onCircle10000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J onCircle10000.txt
	time -o output/Jarvis/J_test_onCircle100000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J onCircle100000.txt
	time -o output/Jarvis/J_test_onCircle1000000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J onCircle1000000.txt
	time -o output/Jarvis/J_test_triangle10.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J triangle10.txt
	time -o output/Jarvis/J_test_triangle100.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J triangle100.txt
	time -o output/Jarvis/J_test_triangle1000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J triangle1000.txt
	time -o output/Jarvis/J_test_triangle10000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J triangle10000.txt
	time -o output/Jarvis/J_test_triangle100000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J triangle100000.txt
	time -o output/Jarvis/J_test_triangle1000000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J triangle1000000.txt

#generates all hull_G_*.txt files with Graham Scan.
.PHONY: runG
runG : convexHull435
	time -o output/Graham/G_test_rectangle10.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G rectangle10.txt
	time -o output/Graham/G_test_rectangle100.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G rectangle100.txt
	time -o output/Graham/G_test_rectangle1000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G rectangle1000.txt
	time -o output/Graham/G_test_rectangle10000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G rectangle10000.txt
	time -o output/Graham/G_test_rectangle100000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G rectangle100000.txt
	time -o output/Graham/G_test_rectangle1000000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G rectangle1000000.txt
	time -o output/Graham/G_test_circle10.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G circle10.txt
	time -o output/Graham/G_test_circle100.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G circle100.txt
	time -o output/Graham/G_test_circle1000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G circle1000.txt
	time -o output/Graham/G_test_circle10000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G circle10000.txt
	time -o output/Graham/G_test_circle100000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G circle100000.txt
	time -o output/Graham/G_test_circle1000000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G circle1000000.txt
	time -o output/Graham/G_test_onCircle10.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G onCircle10.txt
	time -o output/Graham/G_test_onCircle100.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G onCircle100.txt
	time -o output/Graham/G_test_onCircle1000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G onCircle1000.txt
	time -o output/Graham/G_test_onCircle10000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G onCircle10000.txt
	time -o output/Graham/G_test_onCircle100000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G onCircle100000.txt
	time -o output/Graham/G_test_onCircle1000000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G onCircle1000000.txt
	time -o output/Graham/G_test_triangle10.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G triangle10.txt
	time -o output/Graham/G_test_triangle100.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G triangle100.txt
	time -o output/Graham/G_test_triangle1000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G triangle1000.txt
	time -o output/Graham/G_test_triangle10000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G triangle10000.txt
	time -o output/Graham/G_test_triangle100000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G triangle100000.txt
	time -o output/Graham/G_test_triangle1000000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G triangle1000000.txt

#generates all hull_Q_*.txt files with QuickHull.
.PHONY: runQ
runQ : convexHull435
	time -o output/QuickHull/Q_test_rectangle10.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q rectangle10.txt
	time -o output/QuickHull/Q_test_rectangle100.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q rectangle100.txt
	time -o output/QuickHull/Q_test_rectangle1000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q rectangle1000.txt
	time -o output/QuickHull/Q_test_rectangle10000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q rectangle10000.txt
	time -o output/QuickHull/Q_test_rectangle100000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q rectangle100000.txt
	time -o output/QuickHull/Q_test_rectangle1000000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q rectangle1000000.txt
	time -o output/QuickHull/Q_test_circle10.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q circle10.txt
	time -o output/QuickHull/Q_test_circle100.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q circle100.txt
	time -o output/QuickHull/Q_test_circle1000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q circle1000.txt
	time -o output/QuickHull/Q_test_circle10000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q circle10000.txt
	time -o output/QuickHull/Q_test_circle100000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q circle100000.txt
	time -o output/QuickHull/Q_test_circle1000000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q circle1000000.txt
	time -o output/QuickHull/Q_test_onCircle10.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q onCircle10.txt
	time -o output/QuickHull/Q_test_onCircle100.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q onCircle100.txt
	time -o output/QuickHull/Q_test_onCircle1000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q onCircle1000.txt
	time -o output/QuickHull/Q_test_onCircle10000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q onCircle10000.txt
	time -o output/QuickHull/Q_test_onCircle100000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q onCircle100000.txt
	time -o output/QuickHull/Q_test_onCircle1000000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q onCircle1000000.txt
	time -o output/QuickHull/Q_test_triangle10.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q triangle10.txt
	time -o output/QuickHull/Q_test_triangle100.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q triangle100.txt
	time -o output/QuickHull/Q_test_triangle1000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q triangle1000.txt
	time -o output/QuickHull/Q_test_triangle10000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q triangle10000.txt
	time -o output/QuickHull/Q_test_triangle100000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q triangle100000.txt
	time -o output/QuickHull/Q_test_triangle1000000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q triangle1000000.txt

.PHONY: runA
runA : convexHull435
	time -o output/Jarvis/J_test_rectangle10.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J rectangle10.txt
	time -o output/Jarvis/J_test_rectangle100.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J rectangle100.txt
	time -o output/Jarvis/J_test_rectangle1000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J rectangle1000.txt
	time -o output/Jarvis/J_test_rectangle10000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J rectangle10000.txt
	time -o output/Jarvis/J_test_rectangle100000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J rectangle100000.txt
	time -o output/Jarvis/J_test_rectangle1000000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J rectangle1000000.txt
	time -o output/Jarvis/J_test_circle10.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J circle10.txt
	time -o output/Jarvis/J_test_circle100.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J circle100.txt
	time -o output/Jarvis/J_test_circle1000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J circle1000.txt
	time -o output/Jarvis/J_test_circle10000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J circle10000.txt
	time -o output/Jarvis/J_test_circle100000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J circle100000.txt
	time -o output/Jarvis/J_test_circle1000000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J circle1000000.txt
	time -o output/Jarvis/J_test_onCircle10.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J onCircle10.txt
	time -o output/Jarvis/J_test_onCircle100.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J onCircle100.txt
	time -o output/Jarvis/J_test_onCircle1000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J onCircle1000.txt
	time -o output/Jarvis/J_test_onCircle10000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J onCircle10000.txt
	time -o output/Jarvis/J_test_onCircle100000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J onCircle100000.txt
	time -o output/Jarvis/J_test_onCircle1000000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J onCircle1000000.txt
	time -o output/Jarvis/J_test_triangle10.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J triangle10.txt
	time -o output/Jarvis/J_test_triangle100.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J triangle100.txt
	time -o output/Jarvis/J_test_triangle1000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J triangle1000.txt
	time -o output/Jarvis/J_test_triangle10000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J triangle10000.txt
	time -o output/Jarvis/J_test_triangle100000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J triangle100000.txt
	time -o output/Jarvis/J_test_triangle1000000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 J triangle1000000.txt
	time -o output/Graham/G_test_rectangle10.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G rectangle10.txt
	time -o output/Graham/G_test_rectangle100.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G rectangle100.txt
	time -o output/Graham/G_test_rectangle1000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G rectangle1000.txt
	time -o output/Graham/G_test_rectangle10000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G rectangle10000.txt
	time -o output/Graham/G_test_rectangle100000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G rectangle100000.txt
	time -o output/Graham/G_test_rectangle1000000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G rectangle1000000.txt
	time -o output/Graham/G_test_circle10.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G circle10.txt
	time -o output/Graham/G_test_circle100.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G circle100.txt
	time -o output/Graham/G_test_circle1000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G circle1000.txt
	time -o output/Graham/G_test_circle10000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G circle10000.txt
	time -o output/Graham/G_test_circle100000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G circle100000.txt
	time -o output/Graham/G_test_circle1000000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G circle1000000.txt
	time -o output/Graham/G_test_onCircle10.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G onCircle10.txt
	time -o output/Graham/G_test_onCircle100.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G onCircle100.txt
	time -o output/Graham/G_test_onCircle1000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G onCircle1000.txt
	time -o output/Graham/G_test_onCircle10000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G onCircle10000.txt
	time -o output/Graham/G_test_onCircle100000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G onCircle100000.txt
	time -o output/Graham/G_test_onCircle1000000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G onCircle1000000.txt
	time -o output/Graham/G_test_triangle10.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G triangle10.txt
	time -o output/Graham/G_test_triangle100.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G triangle100.txt
	time -o output/Graham/G_test_triangle1000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G triangle1000.txt
	time -o output/Graham/G_test_triangle10000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G triangle10000.txt
	time -o output/Graham/G_test_triangle100000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G triangle100000.txt
	time -o output/Graham/G_test_triangle1000000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 G triangle1000000.txt
	time -o output/QuickHull/Q_test_rectangle10.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q rectangle10.txt
	time -o output/QuickHull/Q_test_rectangle100.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q rectangle100.txt
	time -o output/QuickHull/Q_test_rectangle1000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q rectangle1000.txt
	time -o output/QuickHull/Q_test_rectangle10000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q rectangle10000.txt
	time -o output/QuickHull/Q_test_rectangle100000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q rectangle100000.txt
	time -o output/QuickHull/Q_test_rectangle1000000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q rectangle1000000.txt
	time -o output/QuickHull/Q_test_circle10.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q circle10.txt
	time -o output/QuickHull/Q_test_circle100.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q circle100.txt
	time -o output/QuickHull/Q_test_circle1000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q circle1000.txt
	time -o output/QuickHull/Q_test_circle10000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q circle10000.txt
	time -o output/QuickHull/Q_test_circle100000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q circle100000.txt
	time -o output/QuickHull/Q_test_circle1000000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q circle1000000.txt
	time -o output/QuickHull/Q_test_onCircle10.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q onCircle10.txt
	time -o output/QuickHull/Q_test_onCircle100.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q onCircle100.txt
	time -o output/QuickHull/Q_test_onCircle1000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q onCircle1000.txt
	time -o output/QuickHull/Q_test_onCircle10000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q onCircle10000.txt
	time -o output/QuickHull/Q_test_onCircle100000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q onCircle100000.txt
	time -o output/QuickHull/Q_test_onCircle1000000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q onCircle1000000.txt
	time -o output/QuickHull/Q_test_triangle10.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q triangle10.txt
	time -o output/QuickHull/Q_test_triangle100.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q triangle100.txt
	time -o output/QuickHull/Q_test_triangle1000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q triangle1000.txt
	time -o output/QuickHull/Q_test_triangle10000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q triangle10000.txt
	time -o output/QuickHull/Q_test_triangle100000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q triangle100000.txt
	time -o output/QuickHull/Q_test_triangle1000000.txt -f "Program: %C\nTotal time: %E\n" ./convexHull435 Q triangle1000000.txt
