#include <iostream>
#include <string>
#include <algorithm>
#include <seqan/arg_parse.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/score.h>
#include <seqan/align.h>
#include <seqan/file.h>

using namespace seqan;

struct ModifyStringOptions
{
    CharString inputFile1, inputFile2, primer;
    bool Verbose;
    ModifyStringOptions() :
        Verbose(false)
    {}
};


ArgumentParser::ParseResult
parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("recoverFragment");
    // Set short description, version, and date.
    setShortDescription(parser, "recover fragment into 4 different categories from pair-end seq data");
    setVersion(parser, "0.1");
    setDate(parser, "August 2013");
    
    addOption(parser, seqan::ArgParseOption(
        "I", "inputs", "input of forward and reverse fastq file, path of two files separated by SPACE",
        seqan::ArgParseArgument::STRING, "STR",false, 2));
    setRequired(parser, "I");

    addOption(parser, seqan::ArgParseOption(
        "p", "primer", "fasta file contianing two primer sequences",
        seqan::ArgParseArgument::STRING, "STR"));
    setRequired(parser, "p");

    addOption(parser, seqan::ArgParseOption(
        "v", "verbose", "print alignment information for each alignment"));
   
    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser,
                "\\fBrecoverFragment\\fP \\fB-I\\fP \\fIread_1.fastq\\fP \\fIread_2.fastq\\fP \\fB-p\\fP \\fIprimer.fasta\\fP",
                "store fragment using fasta/fastq into 4 output files 'short_*', 'long_*','evenlong_*','wierd_*'");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
    // Only extract options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;
    // Extract option values.
    getOptionValue(options.inputFile1, parser, "inputs",0);
    getOptionValue(options.inputFile2, parser, "inputs",1);
    getOptionValue(options.primer, parser, "primer",0);
    options.Verbose = isSet(parser,"verbose");

    return seqan::ArgumentParser::PARSE_OK;


}

//local_alignment
template <typename TText>
int local_align(int & score, int & BeginPos1, int & EndPos1, int & BeginPos2, int & EndPos2, TText const & seq1, TText const & seq2, float mismatch, float open_gap, float extend_gap, bool Verbose)
{
    Align<String<char> > ali;
    resize(rows(ali), 2);
    assignSource(row(ali, 0), seq1);
    assignSource(row(ali, 1), seq2);
    score = localAlignment(ali, Score<float, Simple>(2,mismatch,extend_gap, open_gap));
    BeginPos1 = clippedBeginPosition(row(ali, 0));
    EndPos1 = clippedEndPosition(row(ali, 0));
    BeginPos2 = clippedBeginPosition(row(ali, 1));
    EndPos2 = clippedEndPosition(row(ali, 1));
    if (Verbose)
    {
        std::cout << mismatch << "," << open_gap << "," << extend_gap << "\n";
        std::cout << "Score = " << score << std::endl;
        std::cout << ali;
    }
    return 0;

}


int main(int argc, char const ** argv)
{
    // Parse the command line.
    ModifyStringOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    std::cout << "Forward reads file: \t" << options.inputFile1 << '\n'
              << "Reverse reads file: \t" << options.inputFile2 << '\n'
              << "primer fasta file:    \t" << options.primer << '\n';

    // variables
    CharString id,id1,id2,qual1,qual2;
    Dna5String seq1,seq2,frag_seq,primer1,primer2,part_seq1,part_seq2;

    Align<String<char> > ali;
    resize(rows(ali), 2);

    // score1:   seq1         -      seq2_rc
    // score2:   seq1         -      primer1
    // score3:   seq2         -      primer2
    // score4:   trimmed_seq1 - trimmed_seq2
    int score1,score2,score3,score4;
    int BeginPos1, EndPos1, BeginPos2, EndPos2;
    int BeginPos3, EndPos3, BeginPos4, EndPos4;
    int BeginPos5, EndPos5, BeginPos6, EndPos6;
    int n=0, type1_n=0, type2_n=0, type3_n=0, type4_n=0;
    
    // output files
    std::string inputFile1;
    inputFile1 = toCString(options.inputFile1);
    unsigned f1 = inputFile1.find_last_of("/\\");
    if (f1>length(inputFile1)) { f1 = 0; }
    else { f1+=1; }
    unsigned f2 = inputFile1.find_last_of("_");
    std::string type1="short_";
    std::string type2="long_";
    std::string type3_1="evenlong_",type3_2;
    std::string type4_1="wierd_",type4_2;
    type1.append(inputFile1,f1,f2-f1);
    type1.append(".fasta");
    type2.append(inputFile1,f1,f2-f1);
    type2.append(".fasta");
    type3_1.append(inputFile1,f1,f2-f1);
    assign(type3_2,type3_1);
    type3_1.append("_1.fastq");
    type3_2.append("_2.fastq");
    type4_1.append(inputFile1,f1,f2-f1);
    assign(type4_2,type4_1);
    type4_1.append("_1.fastq");
    type4_2.append("_2.fastq");
    SequenceStream Stream_1(type1.c_str(), SequenceStream::WRITE);
    SequenceStream Stream_2(type2.c_str(), SequenceStream::WRITE);
    SequenceStream Stream_3_1(type3_1.c_str(), SequenceStream::WRITE);
    SequenceStream Stream_3_2(type3_2.c_str(), SequenceStream::WRITE);
    SequenceStream Stream_4_1(type4_1.c_str(), SequenceStream::WRITE);
    SequenceStream Stream_4_2(type4_2.c_str(), SequenceStream::WRITE);

    //SequenceStream seqStream()
    // input files
    SequenceStream seqStream1(toCString(options.inputFile1));
    SequenceStream seqStream2(toCString(options.inputFile2));
    SequenceStream seqStream3(toCString(options.primer));
    primer1 = "AGATCGGAAGAGCGGTTCAG";
    primer2 = "AGATCGGAAGAGCGTCGTGT";
    
    readRecord(id, primer1, seqStream3);
    readRecord(id, primer2, seqStream3);

    while (!atEnd(seqStream1))
    {
        /// detect if two fasta has same number of sequences
        if (readRecord(id1, seq1, qual1, seqStream1)!=0 || readRecord(id2, seq2, qual2, seqStream2))
        {
            std::cerr << "ERROR: Could not read from fastq files or they have different numbers of reads !\n";
            return 1;
        }
        n+=1;   
        Dna5String seq2_rc = Dna5StringReverseComplement(seq2);
        local_align(score1, BeginPos1, EndPos1, BeginPos2, EndPos2, seq1, seq2_rc, -2, -4, -1, 
                    options.Verbose);
        if (options.Verbose){
            std::cout << "Aligns Seq1[" << BeginPos1 << ":" << EndPos1 << "]"<<std::endl;
            std::cout << "Aligns Seq2[" << BeginPos2 << ":" << EndPos2 << "]"<<std::endl;
            std::cout << seq1 << std::endl;
            std::cout << seq2_rc << std::endl;
        }
        if (score1 > 1.8*std::max(9,(EndPos1-BeginPos1)) && (BeginPos2 * BeginPos1) == 0)
        {   
            if (EndPos1>=length(seq1) && EndPos1<=length(seq1)+5)
            {
                assign(frag_seq,seq1);
                for (unsigned i = EndPos2-EndPos1+length(seq1); i< length(seq2);++i){
                    appendValue(frag_seq, seq2_rc[i]); }
                writeRecord(Stream_2,id1,frag_seq);
                type2_n+=1;
                if (options.Verbose){
                    std::cout << frag_seq <<std::endl;
                    std::cout << "Type2" << std::endl;}
            } else {
                part_seq1="";
                part_seq2="";
                for (unsigned i = EndPos1; i< length(seq1);++i){
                    appendValue(part_seq1, seq1[i]);
                    appendValue(part_seq2, seq2[i]); }
                local_align(score2, BeginPos3, EndPos3, BeginPos4, EndPos4, part_seq1, primer1, 
                            0,-.5, -.3, options.Verbose);
                local_align(score3, BeginPos5, EndPos5, BeginPos6, EndPos6, part_seq2, primer2, 
                            0,-.5, -.3, options.Verbose);
                if (score2+score3>=1.5*(EndPos5-BeginPos5+EndPos3-BeginPos3)) 
                //use a less stringent cut for the primer match
                {
                    Infix<Dna5String>::Type frag_seq = infix(seq1, BeginPos1, EndPos1);
                    writeRecord(Stream_1,id1,frag_seq);
                    type1_n+=1;
                    if (options.Verbose){
                        std::cout << frag_seq <<std::endl;
                        std::cout << "Type1" << std::endl;}
                } else {
                    writeRecord(Stream_4_1,id1,seq1,qual1);
                    writeRecord(Stream_4_2,id2,seq2,qual2);
                    type4_n+=1;
                    if (options.Verbose){
                        std::cout << "Type4" << std::endl;}
                }    
                
            }
        } else {
            part_seq1="";
            part_seq2="";
            for (unsigned i = 0; i<length(seq2)/3; ++i){
                appendValue(part_seq2, seq2_rc[i]);}
            for (unsigned i = int(2*length(seq1)/3); i<length(seq1); ++i){
                appendValue(part_seq1, seq1[i]); }
            local_align(score4, BeginPos1, EndPos1, BeginPos2, EndPos2, part_seq1, part_seq2,
                            -2, -5, -1, options.Verbose);
            if (score4>=1.7*std::max(7,(EndPos1-BeginPos1)) && BeginPos2==0 && EndPos1>=length(part_seq1)){
                assign(frag_seq,seq1);
                for (unsigned i = EndPos2-BeginPos2; i< length(seq2);++i){
                    appendValue(frag_seq, seq2_rc[i]); }
                writeRecord(Stream_2,id1,frag_seq);
                type2_n+=1;
                if (options.Verbose){
                    std::cout << "Aligns part_Seq2[" << BeginPos2 << ":" << EndPos2 << "]"<<std::endl;
                    std::cout << frag_seq <<std::endl;
                    std::cout << "Type2" << std::endl;}
            } else {
                writeRecord(Stream_3_1,id1,seq1,qual1);
                writeRecord(Stream_3_2,id2,seq2,qual2);
                type3_n+=1;
                if (options.Verbose){
                    std::cout << part_seq1 << std::endl;
                    std::cout << part_seq2 << std::endl;
                    std::cout << "Type3" << std::endl;}
            }
            
        }

    if (n%50000==0)
    {
        std::cout << "Now: " << n << "; Type1: " << type1_n << "; Type2: " << type2_n 
                  << "; Type3: " << type3_n << "; Type4: " << type4_n << ".\n";
    }
    }
    std::cout << "Total: " << n << "; \n";
    std::cout << "\tType1: " << type1_n << "; Type2: " << type2_n << ".\n";
    std::cout << "\tType3: " << type3_n << "; Type4: " << type4_n << ".\n";
    return 0;
}
