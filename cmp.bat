@"C:\Program Files\Java\jdk1.6.0_18\bin\javac" main.java
@pause
@"C:\Program Files\Java\jdk1.6.0_18\bin\jar" cvfm doubleModel.jar mymanifest main.class commonVariables.class commonFunctions.class launcher.class thread.class doubleFundMood.class
@echo -
@del *.class
@pause