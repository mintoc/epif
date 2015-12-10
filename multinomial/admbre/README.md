_multinomialme.tpl_ implements a multinomial mixed effects model in ADMB-RE. ADMB can be obtained [here](http://admb-project.org/).

Having installed ADMB, the code can be compiled with 

**admb -r -s multinomialme**

and ran with 

**./multinomialme -l1 50000000 -l2 200000000 -l3 50000000**

the additional **-l** arguments are to allocate more memory, see [http://bemata.imr.no/admblp/doc/examples/admb_tutorial.html](http://bemata.imr.no/admblp/doc/examples/admb_tutorial.html)
