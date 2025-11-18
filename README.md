# streamlit_demo

Hello! 

This is an example Streamlit app for executing a protein science workflow: designing expression constructs, designing oligonucleotide primers to generate those constructs and generating a picklist to automate dispensing of those primers by a liquid handling robot.

The app is available at [here](https://streamlitdemo-vk4xqig5sq-nw.a.run.app), deployed to Google Cloud Run via GitHub Actions. 

#### Instructions
To use the app:
1. Input the [UniProt](https://www.uniprot.org) ID of a protein of interest. This will load the AlphaFold predicted structure and sequence of the protein of interest.
2. Use the selection slider to choose N an C terminal ends of the protein sequence and click the *add construct boundaries* button to use the selected construct boundaries. The app will generate sequences from every combination of N and C terminal positions entered. The plot at the bottom of the page will show all the constructs generated compared to the full-length protein sequence.
3. Use the sidebar or the *go to primer design page* button to design primers for your protein sequences. A table will display the designed constructs and primers needed to amplify the relevant DNA sequence.
4. Use the sidebar or the *review outputs* button to review the layout of your output plate and download a primer order form (in the format required by the [Merck primer ordering system](https://www.sigmaaldrich.com/GB/en/product/sigma/oligo)) and a picklist to dispense the primers in the format required by an [Echo liquid handler](https://www.mybeckman.uk/landing/ppc/liquid-handlers/acoustic-liquid-handling?utm_source=google&utm_medium=cpc&utm_term=echo+liquid+handler&gad_source=1&gad_campaignid=21016212679&gbraid=0AAAAADsAn-RJVWwVjQGzHFnXyiap_gp0r&gclid=EAIaIQobChMIiJHpsfj7kAMVHohQBh2a9za_EAAYASAAEgLjvPD_BwE).

The app automatically assigns designed constructs a position in a 96-well plate and will warn you if you have designed more constructs than will fit in a plate.

It also uses a simple reverse-translation method to generate an example DNA sequence to use as a primer template. For a production app you would want to use a good codon-optimisation system in place of this step such as the [Twist Biosciences codon optimisation system](https://www.twistbioscience.com/faq/using-your-twist-account/what-does-twist-codon-optimization-tool-do).

