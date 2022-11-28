/*
 *
 * Copyright (C) 2022 Juan Domingo (Juan.Domingo@uv.es)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "jmatrix.h"

// Auxiliary functions:

// Determination of the endianness
unsigned char ThisMachineEndianness()
{
        short int word = 0x0001;
        char *b = (char *)&word;
        return (b[0] ? LITEND : BIGEND);
}

unsigned long long GetFileSize(std::string fname)
{
    	struct stat stat_buf;
    	int rc = stat(fname.c_str(), &stat_buf);
    	if (rc != 0)
    	{
    	 std::string err="Cannot obtain information (with stat system call) of file "+fname+"\n";
    	 err += "This is probably because you are running this in a 32-bit architecture and the file is bigger than 4 GB.\n";
    	 err += "Unfortunately, we have not found yet a solution for that and, if you need to manage so big files,\n";
    	 err += "probably you should consider using a 64-bit architecture.\n";
    	 Rcpp::stop(err);
    	 // NOTE: may be definition of __USE_FILE_OFFSET64 could solve this but it might provoke other problems...
    	}
        return ((unsigned long long)stat_buf.st_size);
}

void PositionsInFile(std::string fname,unsigned long long *start_of_metadata,unsigned long long *start_of_comment)
{
    	struct stat stat_buf;
    	int rc = stat(fname.c_str(), &stat_buf);
    	if (rc != 0)
    	{
    	 std::string err="Cannot obtain information (with stat system call) of file "+fname+"\n";
    	 Rcpp::stop(err);
    	}
    
        *start_of_comment = GetFileSize(fname)-(COMMENT_SIZE+BLOCKSEP_LEN+sizeof(unsigned long long));	
        
    	unsigned long long metadata_pos_mark = GetFileSize(fname)-sizeof(unsigned long long);

        std::ifstream f(fname.c_str());
        f.seekg(metadata_pos_mark,std::ios::beg);
        f.read((char *)start_of_metadata,sizeof(unsigned long long));
        f.close();
        
        return;
}

// Helper functions to generate sensible error messages
std::string MatrixTypeName(unsigned char typeident)
{
    switch (typeident)
    {
        case MTYPEFULL:         return "FullMatrix";
        case MTYPESPARSE:       return "SparseMatrix";
        case MTYPESYMMETRIC:    return "SymmetricMatrix";
        default:                return "UnknownTypeMatrix";
    }
}

std::string DataTypeName(unsigned char dtypeident)
{
    switch (dtypeident)
    {
        case UCTYPE: return "unsigned char";
        case SCTYPE: return "char";
        case USTYPE: return "unsigned short int";
        case SSTYPE: return "short int";
        case UITYPE: return "unsigned int";
        case SITYPE: return "int";
        case ULTYPE: return "unsigned long";
        case SLTYPE: return "long";
        case FTYPE:  return "float";
        case DTYPE:  return "double";
        case LDTYPE: return "long double";
        default: return "Unknown data type";
    }
}


int SizeOfType(unsigned char typeident)
{
    unsigned char vtype=(typeident & 0x0F);
    switch (vtype)
    {
        case NOTYPE:    return -1;
        case UCTYPE:
        case SCTYPE:    return sizeof(char);
        case USTYPE:
        case SSTYPE:    return sizeof(short);
        case UITYPE:
        case SITYPE:    return sizeof(int);
        case ULTYPE:
        case SLTYPE:    return sizeof(long);
        case ULLTYPE:
        case SLLTYPE:   return sizeof(long long);
        case FTYPE:     return sizeof(float);
        case DTYPE:     return sizeof(double);
        case LDTYPE:    return sizeof(long double);
        default:        return -1;
    }
}



