* Encoding: UTF-8.

* Make sure the Active Dataset shows your data file (.sav).
* Steps:
* First. Open HeteroskedasticityV3, select all, and Run Selection (click the triangle green button in the menu above the syntax).
*         You have just invoked the macro in the SPSS memory, which you can call during your SPSS session.
* Second. Highlight the code below, and Run Selection (click the triangle green button above).


BPK dv = q4_sat
       /iv = q8_1_easy q8_2_custom q8_3_deliver 
       /robse = 4.
