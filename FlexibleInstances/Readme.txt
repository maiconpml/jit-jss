Format:
Insput for exetended job shop instances allowing missing operations and recicrulation
Operations are one based, jobs and machines are zero based

O: o[number of operations] J: j[number of jobs] M: m[number of machines]

[processing times, machine zero times should be ommited, order implicit]
sj[size of job in this line] - t i ; [proc.time(index of machine)] [...sj times - one for each opers in job...]
[J lines - one for each job]

Example (job shop with reentrant and zero time operations):

O: 9    J: 3   M: 4

P:
3 - 7 1 ; 5 0 ; 3 1
4 - 1 3 ; 3 1 ; 9 0 ; 3 2
2 - 2 2 ; 1 0

