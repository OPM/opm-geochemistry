/*
 * MIT License
 *
 * Copyright (C) 2025 Aksel Hiorth
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/
#include <opm/simulators/geochemistry/Core/ChemTable.h>

#include <cstddef>

void ChemTable::split_and_add_to_vector(std::vector<std::string> &elements,const std::string &buffer,int skip, const std::string sep)
{
    std::string buf=buffer;
    size_t pos = 0;
    int idx=0;
    while ((pos = buf.find(sep)) != std::string::npos)
    {
        if (idx >= skip)
        {
        elements.push_back(buf.substr(0, pos)); // store the substring
        }
        buf.erase(0, pos + sep.length());  /* erase() function store the current positon and move to next token. */
        idx++;
    }
    if(buf.size()>0) // add last element
        elements.push_back(buf);
}
// give a csv file with a single col value and corresponding row names
// for col_name exists replace values if not create a new col
void ChemTable::add_col_by_name(const std::string col_name, const std::string file_name, const std::string sep)
{
    ChemTable* CM;
    int pos_col = -1;
    CM = new ChemTable;
    CM->read_csv(file_name, sep,false);
    // first find cole name
    auto it = std::find(col_name_.begin(), col_name_.end(), col_name);
    if (it != col_name_.end())
        pos_col = std::distance(col_name_.begin(), it);

    if (pos_col < 0) // column not there
    {
        add_col(std::vector<std::string> {col_name});
        pos_col = noColumns_ - 1;
    }

    int cmi = 0;
    for (auto& row : CM->row_name_)
    {
        auto it2 = std::find(row_name_.begin(), row_name_.end(), row);
        if (it2 != row_name_.end())
        {
            int pos_row = std::distance(row_name_.begin(), it2);
            M_[pos_row][pos_col] = CM->M_[cmi][0];
        }
        cmi++;
    }

    delete CM;

}

void ChemTable::read_csv(const std::string file_name, std::string sep, bool nick_name)
{
    std::vector<std::string>::iterator it;
    // The first entry of database_ is the header
    std::ifstream dbf(file_name, std::ios::in);
    if (!dbf.is_open())
    {
        error_and_exit("ChemTable::read_csv(): Cannot open file \"{:s}\"", file_name);
    }
    std::string buf;
    std::getline(dbf, buf); // read header
    split_and_add_to_vector(col_name_,buf,2,sep);
    noColumns_=col_name_.size();
    int skip_col = 1;
    if (nick_name)
        skip_col = 2;
    while (std::getline(dbf, buf))
    {
        std::vector<std::string> tmp_row;
        split_and_add_to_vector(tmp_row,buf,0,sep);
        it=tmp_row.begin();
        row_name_.push_back(*it);
        std::advance(it,1);
        if (nick_name)
        {
            nick_name_.push_back(*it);
            std::advance(it, 1);
        }
        M_.push_back(std::vector<double>(noColumns_,0.));
        for(;it != tmp_row.end();it++)
        {
            std::size_t idx=it-tmp_row.begin()-skip_col;
            if(idx < M_.back().size())
                M_.back()[idx]=std::stod(*it);
            // else skip elements
        }
    }
    noRows_=M_.size();

}

void ChemTable::write_csv(const std::string file_name, [[maybe_unused]] std::string sep)
{
    std::ofstream dbf(file_name, std::ios::out);
    if (!dbf.is_open())
    {
        error_and_exit("ChemTable::write_csv(): Cannot write to file \"{:s}\"", file_name);
    }

    dbf <<"Name\tNickName" <<std::setprecision(output_precision_);
    for(int i = 0 ;i< noColumns_;++i)
        dbf<< "\t"<<col_name_[i];
    dbf <<"\n";
    for (int j = 0;j < noRows_; ++j)
        {
            dbf <<row_name_[j]<<"\t"<<nick_name_[j];
            for(int i = 0 ;i< noColumns_;++i)
            {
                dbf << "\t"<<M_[j][i];
            }
            dbf <<"\n";
        }
    dbf.close();
}
// drops columns in M_ - matrix and col_name_ that are listed in std::vector
void ChemTable::drop_col(std::vector<int> &columns_to_drop)
{
    if (columns_to_drop.size() != static_cast<std::size_t>(noColumns_))
        error_and_exit("ChemTable::drop_col(): Enter one value for each column {:d}", columns_to_drop.size());
    std::vector<std::vector<double>> M_new;
    std::vector<std::string> col_name_new;
    for(int i=0;i<noColumns_;++i)
        if(!columns_to_drop[i])
            col_name_new.push_back(col_name_[i]);
    for (int j = 0;j < noRows_; ++j)
        {
       M_new.push_back(std::vector<double>());
            for(int i=0;i<noColumns_;++i)
                {
                 if(!columns_to_drop[i])
                     M_new.back().push_back(M_[j][i]);
                }
        }
    noColumns_=col_name_new.size();
    M_=M_new;
    col_name_=col_name_new;
}
// append only columns with exact same name
void ChemTable::append(const ChemTable& CM)
{
    std::vector<int> map_col;
    bool found=false;

    for (int j = 0; j < noColumns_; ++j)
    {
        map_col.push_back(-1);
        for (int i =0;i<CM.noColumns_;++i)
            if(CM.col_name_[i]==col_name_[j])
                {
                map_col[j]=i;
                found=true;
                }
    }
    if (!found) // no common coloumns
        return;

    for (int j = 0; j < CM.noRows_; ++j)
    {
        row_name_.push_back(CM.row_name_[j]);
        nick_name_.push_back(CM.nick_name_[j]);
        M_.push_back(std::vector<double>(noColumns_,0.));
        for (int i = 0; i < noColumns_; ++i)
        {
            if (map_col[i] > -1)
            {
                M_.back()[i]=CM.M_[j][map_col[i]];
            }
        }
    }
    noRows_=row_name_.size();
}
/* renames first occourance, if col not found - do nothing*/
void ChemTable::rename_col(std::string old_name, std::string new_name)
{

    for (auto it = col_name_.begin(); it != col_name_.end(); ++it)
    {
        if (*it == old_name)
        {
            *it = new_name;
            break;
        }
    }
}

/* renames first occourance, if col not found - do nothing*/
void ChemTable::rename_col(const std::vector<std::string> &old_names, const std::vector<std::string> &new_names)
{
    if (old_names.size() != new_names.size())
        error_and_exit("ChemTable::rename_col(): Enter one old and one new name");
    for (auto it = old_names.begin(); it != old_names.end(); ++it)
    {
        int idx = it - old_names.begin();
        rename_col(old_names[idx], new_names[idx]);
    }
}

void ChemTable::swap_col(std::string col_name, int new_pos)
{
    auto it = find(col_name_.begin(), col_name_.end(), col_name);
    if (it != col_name_.end())
    {
        int idx = it - col_name_.begin();
        double values_in_old_place;
        std::string col_name_in_old_place = col_name_[new_pos];
        col_name_[new_pos] = *it;
        col_name_[idx] = col_name_in_old_place;
        for (int i = 0; i < noRows_; ++i)
        {
            values_in_old_place = M_[i][new_pos]; // save old value;
            M_[i][new_pos] = M_[i][idx];
            M_[i][idx] = values_in_old_place;
        }

    }
}

void ChemTable::move_col(std::string col_name, int new_pos)
{
    auto it = find(col_name_.begin(), col_name_.end(), col_name);
    if (it != col_name_.end())
    {
        int idx = it - col_name_.begin();


        if (idx == new_pos) return;
        if (idx < new_pos)
        {
            std::string tmp = col_name_[idx];
            for (int j = idx; j < new_pos; ++j)
                col_name_[j] = col_name_[j + 1];
            col_name_[new_pos] = tmp;
            for (int i = 0; i < noRows_; ++i)
            {
                double values_in_old_place = M_[i][idx];
                for (int j = idx; j < new_pos; ++j)
                {
                    M_[i][j] = M_[i][j + 1];
                }
                M_[i][new_pos] = values_in_old_place;
            }
        }
        else
        {
            std::string tmp = col_name_[idx];
            for (int j = idx; j > new_pos; --j)
                col_name_[j] = col_name_[j - 1];
            col_name_[new_pos] = tmp;
            for (int i = 0; i < noRows_; ++i)
            {
                double values_in_old_place = M_[i][idx];
                for (int j = idx; j > new_pos; --j)
                {
                    M_[i][j] = M_[i][j - 1];
                }
                M_[i][new_pos] = values_in_old_place;
            }
        }
    }
}
void ChemTable::sort_alphabetically(bool reverse)
{
    std::vector<std::string> row_name_new, nick_name_new;
    std::vector<std::vector<double>> Mnew;
    row_name_new.resize(row_name_.size());
    nick_name_new.resize(nick_name_.size());
    Mnew.resize(M_.size());
    int idx = 0;
    for (int i : sort_vector_indices(row_name_,reverse))
    {
        row_name_new[idx] = row_name_[i];
        nick_name_new[idx] = nick_name_[i];
        Mnew[idx] = M_[i];
        ++idx;
    }
    row_name_ = row_name_new;
    nick_name_ = nick_name_new;
    M_ = Mnew;
}
std::unique_ptr<ChemTable> ChemTable::create_reduced_table(
    const std::vector<int>& columns_to_keep,
    int start_row, int start_col
) const
{
    auto newChemTable = std::make_unique<ChemTable>();
    ChemTable& CM = *newChemTable;

    int N_col = noColumns_ - start_col + 1;
    for(int i=0; i < noColumns_ - start_col + 1; ++i)
    {
        if(columns_to_keep[i] == 0)     --N_col;
    }

    // Remove rows that are not needed
    int N_row = noRows_ - start_row + 1;
    for(int i=0; i < noRows_ - start_row + 1; ++i)
    {
        for(int j=0; j < noColumns_ - start_col + 1; ++j)
        {
            const bool row_depends_on_column = ( M_[start_row + i - 1][start_col + j - 1] != 0.0 );
            if(columns_to_keep[j] == 0 && row_depends_on_column)
            {
                --N_row;
                break;
            }
        }
    }

    CM.noRows_ = N_row;
    CM.noColumns_ = N_col;

    // Allocate sufficient space
    CM.row_name_.resize(N_row);
    CM.nick_name_.resize(N_row);
    CM.col_name_.resize(N_col);
    CM.abs_pos_.resize(N_row);
    CM.abs_pos_bas_.resize(N_col);

    // Allocate space for matrix
    CM.M_.resize(CM.noRows_);
    for (int i = 0; i < CM.noRows_; ++i)
    {
        CM.M_[i].resize(CM.noColumns_);
    }

    // Concentration of basis species
    CM.log_m_.resize(CM.noRows_);
    CM.log_a_.resize(CM.noRows_);
    CM.log_g_.resize(CM.noRows_);
    CM.log_QK_.resize(CM.noRows_);
    CM.delta_.resize(CM.noRows_);

    // Initialize concentrations
    for(int i=0; i < CM.noRows_; ++i)
    {
        CM.log_QK_[i] = 0.0;
        CM.delta_[i] = 0.0;
        CM.log_g_[i] = 0.0;
        CM.log_m_[i] = NumericalConstants::LOG10_ALMOST_ZERO;
        CM.log_a_[i] = NumericalConstants::LOG10_ALMOST_ZERO;
    }

    int Npos_j = 0;
    for(int j=0; j < noColumns_ - start_col + 1 ; ++j)
    {
        if(columns_to_keep[j] == 1)
        {
            CM.col_name_[Npos_j] = col_name_[start_col + j - 1];
            CM.abs_pos_bas_[Npos_j] = j;
            ++Npos_j;
        }
    }

    int Npos_i = 0;

    for(int i=0; i < noRows_ - start_row + 1; ++i)
    {
        Npos_j = 0;
        bool rowIsIncluded = true;

        for(int j=0; j < noColumns_ - start_col + 1 ; ++j)
        {
            const bool row_depends_on_column = ( M_[start_row + i - 1][start_col + j - 1] != 0.0 );
            if(columns_to_keep[j] == 0 && row_depends_on_column)
            {
                --Npos_i;
                rowIsIncluded = false;
                break;
            }
        }

        if(rowIsIncluded)
        {
            for(int j=0; j < noColumns_ - start_col + 1 ; ++j)
            {
                if(columns_to_keep[j] == 1)
                {
                    CM.abs_pos_[Npos_i] = i;
                    CM.M_[Npos_i][Npos_j] = M_[start_row + i - 1][start_col + j - 1];
                    CM.row_name_[Npos_i] = row_name_[start_row + i - 1];
                    CM.nick_name_[Npos_i] = nick_name_[start_row + i - 1];
                    ++Npos_j;
                }
            }
        }
        ++Npos_i;
    }
    copy_data_to_new_table(CM);

    CM.logK_.resize(CM.noRows_, 0.0);
    CM.dG_TP_.resize(CM.noRows_, 0.0);
    CM.re_.resize(CM.noRows_, 0.0);

    CM.update_number_of_analytical_species();
    CM.set_up_sparse_matrix();

    return newChemTable;
}

std::unique_ptr<ChemTable> ChemTable::create_reduced_table_rc(
    const std::vector<int>& rows_to_keep,
    const std::vector<int>& columns_to_keep,
    int start_row,
    int start_col
) const
{
    auto newChemTable = std::make_unique<ChemTable>();
    ChemTable& CM = *newChemTable;

    int N_col = noColumns_ - start_col + 1;
    for(int i=0; i < noColumns_ - start_col + 1; ++i)
    {
        if(columns_to_keep[i] == 0)     --N_col;
    }

    // Remove rows that are not needed
    int N_row = noRows_ - start_row + 1;
    for(int i=0; i < noRows_ - start_row + 1; ++i)
    {
        if(rows_to_keep[i] == 0)     --N_row;
    }

    if(N_row == 0 && N_col == 0)    return newChemTable;  // Not really needed..?

    CM.noRows_ = N_row;
    CM.noColumns_ = N_col;

    // Allocate sufficient space
    CM.row_name_.resize(N_row);
    CM.nick_name_.resize(N_row);
    CM.col_name_.resize(N_col);
    CM.abs_pos_.resize(N_row);
    CM.abs_pos_bas_.resize(N_col);

    // Allocate space for matrix
    CM.M_.resize(CM.noRows_);
    for (int i = 0; i < CM.noRows_; ++i)
    {
        CM.M_[i].resize(CM.noColumns_);
    }

    // Concentration of basis species
    CM.log_m_.resize(CM.noRows_);
    CM.log_a_.resize(CM.noRows_);
    CM.log_g_.resize(CM.noRows_);
    CM.log_QK_.resize(CM.noRows_);
    CM.delta_.resize(CM.noRows_);

    // Initialize concentrations
    for(int i=0; i < CM.noRows_; ++i)
    {
        CM.log_QK_[i] = 0.0;
        CM.delta_[i]  = 0.0;
        CM.log_g_[i]  = 0.0;
        CM.log_m_[i]  = NumericalConstants::LOG10_ALMOST_ZERO;
        CM.log_a_[i]  = NumericalConstants::LOG10_ALMOST_ZERO;
    }

    // 13/9-22: Set to zero here, otherwise it is always zero when used inside
    //          the loop.
    int Npos_j = 0;
    for(int j=0; j < noColumns_ - start_col + 1 ; ++j)
    {
        if(columns_to_keep[j] == 1)
        {
            CM.col_name_[Npos_j] = col_name_[start_col + j - 1];
            CM.abs_pos_bas_[Npos_j] = j;
            ++Npos_j;
        }
    }

    int Npos_i = 0;
    for(int i=0; i < noRows_ - start_row + 1; ++i)
    {
        Npos_j = 0;
        bool add_column_and_row_element = false;

        for(int j=0; j < noColumns_ - start_col + 1 ; ++j)
        {
            if( (columns_to_keep[j] == 1) && (rows_to_keep[i] == 1 ) )
            {
                add_column_and_row_element = true;
                CM.abs_pos_[Npos_i] = i;
                CM.M_[Npos_i][Npos_j] = M_[start_row + i - 1][start_col + j - 1];
                CM.row_name_[Npos_i] = row_name_[start_row + i - 1];
                CM.nick_name_[Npos_i] = nick_name_[start_row + i - 1];
                ++Npos_j;
            }
        }

        if(add_column_and_row_element)    ++Npos_i;
    }
    // There are a few vectors that we do not actually need for the case of
    // the basis species. However, we still allocate memory for them, as this
    // lets us reuse the same code for basis species, secondary complexes,
    // and minerals.
    copy_data_to_new_table(CM);

    CM.logK_.resize(CM.noRows_, 0.0);
    CM.dG_TP_.resize(CM.noRows_, 0.0);
    CM.re_.resize(CM.noRows_, 0.0);

    CM.update_number_of_analytical_species();
    CM.set_up_sparse_matrix();

    return newChemTable;
}

void ChemTable::write(const std::string& file_name) const
{
    std::ofstream fpd(file_name, std::ios::binary);
    auto buffer_row = fmt::memory_buffer();

    // Write header -->
    format_to(std::back_inserter(buffer_row), "row_name\t");
    if (!model_.empty()) format_to(std::back_inserter(buffer_row), "model\t");
    if (!type_.empty()) format_to(std::back_inserter(buffer_row), "type\t");
    if (!logK_.empty()) format_to(std::back_inserter(buffer_row), "logK\t");
    if (!mol_volume_.empty()) format_to(std::back_inserter(buffer_row), "mol_volume\t");
    if (!mol_weight_.empty()) format_to(std::back_inserter(buffer_row), "mol_weight\t");
    if (!log_a_.empty()) format_to(std::back_inserter(buffer_row), "log_a\t");
    if (!log_m_.empty()) format_to(std::back_inserter(buffer_row), "log_m\t");
    if (!log_g_.empty()) format_to(std::back_inserter(buffer_row), "log_g\t");
    if (!log_af_.empty()) format_to(std::back_inserter(buffer_row), "log_af\t");
    if (!fugacity_.empty()) format_to(std::back_inserter(buffer_row), "fugacity\t");
    if (!delta_.empty()) format_to(std::back_inserter(buffer_row), "delta\t");
    if (!log_QK_.empty()) format_to(std::back_inserter(buffer_row), "log_QK\t");
    if (!a0_.empty()) format_to(std::back_inserter(buffer_row), "a0\t");
    if (!charge_.empty()) format_to(std::back_inserter(buffer_row), "charge\t");
    if (!scharge_.empty()) format_to(std::back_inserter(buffer_row), "scharge\t");
    if (!deltaG_.empty()) format_to(std::back_inserter(buffer_row), "DeltaG\t");
    if (!deltaH_.empty()) format_to(std::back_inserter(buffer_row), "DeltaH\t");
    if (!S_.empty()) format_to(std::back_inserter(buffer_row), "S\t");
    if (!a1_.empty()) format_to(std::back_inserter(buffer_row), "a1\t");
    if (!a2_.empty()) format_to(std::back_inserter(buffer_row), "a2\t");
    if (!a3_.empty()) format_to(std::back_inserter(buffer_row), "a3\t");
    if (!a4_.empty()) format_to(std::back_inserter(buffer_row), "a4\t");
    if (!c1_.empty()) format_to(std::back_inserter(buffer_row), "c1\t");
    if (!c2_.empty()) format_to(std::back_inserter(buffer_row), "c2\t");
    if (!dG_TP_.empty()) format_to(std::back_inserter(buffer_row), "deltaG_TP\t");

    for(const auto& col_name: col_name_) format_to(std::back_inserter(buffer_row), "{:s}\t", col_name);
    format_to(std::back_inserter(buffer_row), "\n");

    auto header = std::string(buffer_row.begin(), buffer_row.end());
    fmt::print(fpd, "{}", header);
    // <--- Done with header

    for(int i=0; i < noRows_; ++i)
    {
        buffer_row.clear();

        format_to(std::back_inserter(buffer_row), "{:s}\t", row_name_[i]);
        if (!model_.empty()) format_to(std::back_inserter(buffer_row), "{:d}\t", static_cast<int>(model_[i]));
        if (!type_.empty()) format_to(std::back_inserter(buffer_row), "{:d}\t", type_[i]);
        if (!logK_.empty()) format_to(std::back_inserter(buffer_row), "{:f}\t", logK_[i]);
        if (!mol_volume_.empty()) format_to(std::back_inserter(buffer_row), "{:f}\t", mol_volume_[i]);
        if (!mol_weight_.empty()) format_to(std::back_inserter(buffer_row), "{:f}\t", mol_weight_[i]);
        if (!log_a_.empty()) format_to(std::back_inserter(buffer_row), "{:f}\t", log_a_[i]);
        if (!log_m_.empty()) format_to(std::back_inserter(buffer_row), "{:f}\t", log_m_[i]);
        if (!log_g_.empty()) format_to(std::back_inserter(buffer_row), "{:f}\t", log_g_[i]);
        if (!log_af_.empty()) format_to(std::back_inserter(buffer_row), "{:f}\t", log_af_[i]);
        if (!fugacity_.empty()) format_to(std::back_inserter(buffer_row), "{:f}\t", fugacity_[i]);
        if (!delta_.empty()) format_to(std::back_inserter(buffer_row), "{:f}\t", delta_[i]);
        if (!log_QK_.empty()) format_to(std::back_inserter(buffer_row), "{:f}\t", log_QK_[i]);
        if (!a0_.empty()) format_to(std::back_inserter(buffer_row), "{:f}\t", a0_[i]);
        if (!charge_.empty()) format_to(std::back_inserter(buffer_row), "{:f}\t", charge_[i]);
        if (!scharge_.empty()) format_to(std::back_inserter(buffer_row), "{:f}\t", scharge_[i]);
        if (!deltaG_.empty()) format_to(std::back_inserter(buffer_row), "{:f}\t", deltaG_[i]);
        if (!deltaH_.empty()) format_to(std::back_inserter(buffer_row), "{:f}\t", deltaH_[i]);
        if (!S_.empty()) format_to(std::back_inserter(buffer_row), "{:f}\t", S_[i]);
        if (!a1_.empty()) format_to(std::back_inserter(buffer_row), "{:f}\t", a1_[i]);
        if (!a2_.empty()) format_to(std::back_inserter(buffer_row), "{:f}\t", a2_[i]);
        if (!a3_.empty()) format_to(std::back_inserter(buffer_row), "{:f}\t", a3_[i]);
        if (!a4_.empty()) format_to(std::back_inserter(buffer_row), "{:f}\t", a4_[i]);
        if (!c1_.empty()) format_to(std::back_inserter(buffer_row), "{:f}\t", c1_[i]);
        if (!c2_.empty()) format_to(std::back_inserter(buffer_row), "{:f}\t", c2_[i]);
        if (!dG_TP_.empty()) format_to(std::back_inserter(buffer_row), "{:f}\t", dG_TP_[i]);

        if(!M_.empty())
        {
            for(int j=0; j < noColumns_; ++j)
            {
                format_to(std::back_inserter(buffer_row), "{:f}\t", M_[i][j]);
            }
        }
        format_to(std::back_inserter(buffer_row), "\n");
        auto row_str = std::string(buffer_row.begin(), buffer_row.end());
        fmt::print(fpd, "{}", row_str);
    }

    if(!T_.empty())
    {
        fmt::print(fpd, "T\n");
        for(const auto& temperature : T_) fmt::print(fpd, "{:f}\t", temperature);
        fmt::print(fpd, "\n");
    }
}

/*
 * Searches for name in row - returns position index where it is found, or -1 if non-existent.
 * Is not case-sensitive.
 *
 * NB: For species other than E and e-, also searches for the nicknames!
 */
int ChemTable::get_row_index(const std::string& name) const
{
    const std::string rw = to_lower_case(name);

    for (int j = 0; j < noRows_; ++j)
    {
        const std::string rw_j = to_lower_case(row_name_[j]);
        const std::string nw_j = to_lower_case(nick_name_[j]);
        // HACK (23/2-22): The nicknames of E and e- are the same!
        //                 At least for now, manually prevent this comparison to be true...
        const bool ad_hoc_nick_name_comparison = (rw != "e") && (nw_j == rw);
        if (rw_j == rw || ad_hoc_nick_name_comparison)
        {
            return j;
        }
    }
    return -1;
}

/*
 * Searches for name in column - returns position index where it is found, or -1 if non-existent.
 * Is not case-sensitive.
 */
int ChemTable::get_col_index(const std::string& name) const
{
    const std::string cw = to_lower_case(name);

    for (int j = 0; j < noColumns_; ++j)
    {
        std::string cw_j = to_lower_case(col_name_[j]);
        if (cw_j == cw )
        {
            return j;
        }
    }
    return -1;
}

/* Returns true if the input name, or the associated nickname, is in found
 * in the given row. Is not case-sensitive. */
bool ChemTable::name_is_in_row(const std::string& name, int row_index) const
{
    const std::string rw = to_lower_case(name);

    const std::string rw_j = to_lower_case(row_name_[row_index]);
    const std::string nw_j = to_lower_case(nick_name_[row_index]);

    if (rw_j == rw || nw_j == rw) return true;
    else return false;
}

/**
 * @return A vector with values 1 for elements that exist and 0 otherwise. */
std::vector<int> ChemTable::find_row_elements(const std::vector<std::string>& names) const
{
    std::vector<int> val;
    val.resize(row_name_.size(), 0);

    for (auto it = row_name_.cbegin(); it != row_name_.cend(); ++it)
    {
        const auto pos = static_cast<int>(it - row_name_.cbegin());

        for(const auto& name: names)
        {
            std::string NAME = to_upper_case(name);
            if (NAME == to_upper_case(*it) || NAME == to_upper_case(nick_name_[pos]))
            {
                val[pos] = 1;
                break;
            }
        }
    }

    return val;
}

/** @return The number of duplicate entries. */
int ChemTable::check_for_duplicate_species(bool is_case_sensitive, bool check_nick_names)
{
    const std::string start_msg = check_nick_names ? "Duplicate nick name for" : "Duplicate";

    auto names_to_check = check_nick_names ? nick_name_ : row_name_;
    if(!is_case_sensitive) names_to_check = to_lower_case(names_to_check);

    const auto& duplicate_species = find_duplicates(names_to_check);
    for(const auto& species: duplicate_species){
        std::cout << start_msg << " species: " << species << "\n";
    }
    return static_cast<int>(duplicate_species.size());
}

bool ChemTable::table_has_duplicate_row_names(bool case_sensitive) const{

    if(case_sensitive){
        return vector_has_duplicates(row_name_);
    }
    return vector_has_duplicates(to_lower_case(row_name_));
}

bool ChemTable::table_has_duplicate_nick_names(bool case_sensitive) const{

    if(case_sensitive){
        return vector_has_duplicates(nick_name_);
    }
    return vector_has_duplicates(to_lower_case(nick_name_));

}

void ChemTable::initialize_full_database_hkf(ChemTable::Type type, const std::vector<std::pair<std::string,std::string>> &replace_names )
{

    if(type == ChemTable::Type::Basis)
    {
        database_ = create_basis_species_table_hkf();
        add_surface_species_to_basis_table(database_);
    }
    else if(type == ChemTable::Type::Complex)
    {
        database_ = create_complex_species_table_hkf();
        add_surface_species_to_complex_table(database_);
    }
    else if(type == ChemTable::Type::Mineral)
    {
        database_ = create_minerals_table_hkf();
        add_surface_species_to_minerals_table(database_);
    }
    else
    {
        error_and_exit("ChemTable must hold either basis species, secondary species, or buffers.");
    }

    // The first entry of database_ is the header
    noRows_ = database_.size() - 1;

    std::stringstream buf;
    buf.str(database_[0]);

    std::string bn;
    buf >> bn; // Name
    buf >> bn; // NickName
    while (buf >> bn)
    {
        col_name_.push_back(bn);
    }

    noColumns_ = col_name_.size();

    M_.resize(noRows_);
    for (std::size_t i = 0; i < database_.size()-1; ++i)
    {
        buf.clear();
        buf.str(database_[i+1]);
        buf >> bn;
        row_name_.push_back(bn);
        buf >> bn;
        nick_name_.push_back(bn);

        M_[i].resize(noColumns_, 0.0);

        double mi;
        int count = 0;
        while (buf >> mi)
        {
            M_[i][count] = mi;
            ++count;
        }

        if (count != noColumns_)
        {
            std::cout << "Could not parse all floating-point values for type=\"" << typeToString(type) << "\", you should check the database definition!\n";
        }
    }

    model_.resize(noRows_, LogKModel::HKF);
    for (int i = 0; i < size_analytical_; ++i)
    {
        model_[noRows_ - size_analytical_] = LogKModel::ANA;  // place at the end
    }

    // Change some of the Nick names
    if (!replace_names.empty())
    {
        for (auto it = replace_names.begin(); it != replace_names.end(); ++it)
        {
            int pos = in_list((*it).first, nick_name_);
            if (pos > -1)
            {
                nick_name_[pos] = (*it).second;
            }
        }
    }
}

void ChemTable::init_tables()
{
    type_       = pick_column<int>("type");
    logK_       = pick_column<double>("logK");
    log_af_     = pick_column<double>("logaf");
    a0_         = pick_column<double>("a0");
    charge_     = pick_column<double>("charge");
    scharge_    = pick_column<double>("scharge");
    mol_volume_ = pick_column<double>("mol_volume");
    mol_weight_ = pick_column<double>("mol_weight");

    //  Name    NickName    DeltaG  DeltaH  S   a1  a2  a3  a4  c1  c2  omega
    deltaG_     = pick_column<double>("DeltaG");
    deltaH_     = pick_column<double>("DeltaH");
    S_          = pick_column<double>("S");
    a1_         = pick_column<double>("a1");
    a2_         = pick_column<double>("a2");
    a3_         = pick_column<double>("a3");
    a4_         = pick_column<double>("a4");
    c1_         = pick_column<double>("c1");
    c2_         = pick_column<double>("c2");
    a7_         = pick_column<double>("a7");
    omega_      = pick_column<double>("omega");

    // These vectors hold HKF properties derived later from the database rows.
    re_.resize(noRows_, 0.0);
    dG_TP_.resize(noRows_, 0.0);
    molar_volume_PT_.resize(noRows_, 0.0);
    Cp_PT_.resize(noRows_, 0.0);
}

void  ChemTable::add_cubic_EOS_parameters_for_gases()
{
    std::vector<std::string> eos_data;
    eos_data = add_EOS_data_gases();
    std::string name;
    double Tc;
    double Pc;
    double Omega;
    omega_acc_.resize(noRows_, 0.);
    Pc_.resize(noRows_, 0.);
    Tc_.resize(noRows_, 0.);

    for (const auto& line : eos_data)
    {

        std::stringstream buf(line);
        buf >> name >> Tc >> Pc >> Omega;
        int pos=get_row_index(name);
        if (pos > 0)
        {
            Pc_[pos] = Pc;
            Tc_[pos] = Tc;
            omega_acc_[pos] = Omega;
        }

    }
}


void ChemTable::remove_row_species(const std::vector<std::string>& species_to_remove)
{
    if(species_to_remove.empty()) return;
    else erase_row_elements(species_to_remove);
}

/** Adds new  column species to the database, hence (potentially) resizes the number of columns.
 *  Returns the number of species added.
 */
int ChemTable::add_col(const std::vector<std::string>& names, std::vector<double> default_values)
{
    if (!default_values.empty())
    {
        if (names.size() != default_values.size())
        {
            error_and_exit("add_col: default_values has to match column names!");
        }
    }
    int no_species_added = 0;
    double default_d = 0.;
    int i = 0;
    for(const auto& name: names)
    {
        if (!default_values.empty())
        {
            default_d = default_values[i];
            ++i;
        }

        if(!element_is_in_container(name, col_name_))
        {
            // Add a column for the new species
            resize(noRows_, noColumns_ + 1,default_d);
            col_name_[noColumns_ - 1].assign(name);
            ++no_species_added;
        }
    }
    return no_species_added;
}

/** Adds new  row species to the database, hence (potentially) resizes the number of rows..
 *  Returns the number of species added.
 */
int ChemTable::add_species_row(const std::string& name, const std::vector<double>& input_values)
{
    int no_added_species = 0;

    int row_index = -1;
    bool found_specie = false;
    for (int i = 0; i < noRows_; ++i)
    {
        found_specie = name_is_in_row(name, i);
        if (found_specie)
        {
            row_index = i;
            ++no_added_species;
            break;
        }
    }
    if(!found_specie) assert(row_index == -1);

    // If element already exists in db, erase the old one (more flexible)
    if (found_specie) erase_row(row_index);

    assert(noRows_ >= 1);  // Not in db (?)

    resize(noRows_ + 1, noColumns_);
    row_name_[noRows_ - 1] = name;
    nick_name_[noRows_ - 1] = name;

    // i) Fill a temporary vector with the input values. If any values are
    //    missing in the input given by the user, repeat the last number
    //    as many times as needed.
    //
    // ii) Place the new values in the last row. This is ok, because we
    //     just made room for one more.
    //
    if (!input_values.empty())  // only relevant for basis species...
    {
        std::vector<double> temp;
        temp.resize(noColumns_, input_values.back());

        std::size_t no_values = noColumns_;
        if (input_values.size() < no_values)    no_values = input_values.size();

        for (std::size_t i = 0; i < no_values; ++i)     temp[i] = input_values[i];

        for (int i = 0; i < noColumns_; ++i)    M_[noRows_ - 1][i] = temp[i];
    }

    return no_added_species;  // Can only be 0 or 1...
}

/** Add stoichiometric reaction. */
int ChemTable::add_reaction(const std::string& reac,
                            int type,
                            LogKModel model,
                            double charge,
                            double a0,
                            double mol_volume,
                            const std::vector<double>& hkf,
                            const std::vector<std::string>& basis_specie_names,
                            const std::vector<double>& stoic)
{
    std::vector<int> pos;
    for (const auto& bas_name : basis_specie_names)
    {
        const int idx = index_of_element_in_container(bas_name, col_name_);
        if (idx >= 0) pos.push_back(idx);  // valid basis species
        else{
            std::cout << "Error in geochem::add_species: I am not adding the following reaction: ";
            std::cout << reac << ", because basis species " << bas_name << " does not exist\n";
            return -1;
        }
    }

    const int j = index_of_element_in_container(reac, row_name_);
    if(j == -1)
    {
        std::cout << "Could not find element " << reac << ", do not add it to the database \n";
        return -1;
    }

    // If we get to this point, we know the element exists in db
    model_[j] = model;
    for (std::size_t i = 0; i < stoic.size(); ++i)
    {
        M_[j][pos[i]] = stoic[i];
    }

    int k = index_of_element_in_container("type", col_name_);
    if(k >= 0) M_[j][k] = type;

    k = index_of_element_in_container("a0", col_name_);
    if(k >= 0)  M_[j][k] = a0;

    k = index_of_element_in_container("charge", col_name_);
    if(k >= 0 && (type == GeochemicalComponentType::AQUEOUS_COMPLEX || type == GeochemicalComponentType::ION_EXCHANGE))
    {
        M_[j][k] = charge;
    }
    else if(k >= 0) M_[j][k] = 0.0;

    k = index_of_element_in_container("scharge", col_name_);
    if(k >= 0)
    {
        M_[j][k] = (type == GeochemicalComponentType::SURFACE_COMPLEX) ? charge : 0.0;
    }

    k = index_of_element_in_container("mol_volume", col_name_);
    if(k >= 0) M_[j][k] = mol_volume;

    k = index_of_element_in_container("DeltaG", col_name_);
    if(k >= 0)
    {
        for (std::size_t offset = 0; offset < hkf.size(); ++offset)
        {
            M_[j][k+offset] = hkf[offset];  // DeltaG, ... (and the rest)
        }
    }

    return 1;
}

/**
 * @return The number of basis species in mineral, not counting elements in the ignore list.
 */
int ChemTable::find_no_basis_species(int pos_buf, int* pos_bas, const std::vector<int>& ignore_list)
{
    int species_count = 0;
    for (int i = noColumns_ - 1; i >= 0; --i)
    {
        if(M_[pos_buf][i] != 0 && !element_is_in_container(i, ignore_list))
        {
            ++species_count;
            *pos_bas = i;  // Can be updated multiple times, if so we use the final value (== the lowest column index).
        }
    }
    return species_count;
}

/*
 in_fmatrix:
 from a matrix extract the rows and columns given in the list
 INPUT: 1) fmatrix: input matrix to be shaved
 2) row: vector containing a list of integers, row number 2, 4, etc.
 the integers don't have to be sorted
 3) dim_row: dimension of the row vector
 4) col: vector containing a list of integers, col number 1, 2, etc.
 5) dim_col: dimension of the col vector
 OUTPUT: 1) f_red: the reduced matrix
 RETURN: void
 */
void ChemTable::in_fmatrix(const int* row, int dim_row, const int* col, int dim_col, double** f_red)
{
    for(int i=0; i < dim_row; ++i)
    {
        const int row_element = row[i];
        for(int j=0; j < dim_col; ++j)
        {
            const int col_element = col[j];
            f_red[i][j] = M_[row_element][col_element];
        }
    }
}

void ChemTable::resize(int no_rows, int no_columns,double default_col_val)
{
    if (no_rows != noRows_)         resize_row(no_rows);
    if (no_columns != noColumns_)   resize_col(no_columns,default_col_val);
}

void ChemTable::transpose()
{
    M_trans_.resize(noColumns_);
    for (int i = 0; i < noColumns_; ++i)
    {
        M_trans_[i].resize(noRows_);
    }

    for (int i = 0; i < noRows_; ++i)
    {
        for (int j = 0; j < noColumns_; ++j)
        {
            M_trans_[j][i] = M_[i][j];
        }
    }
}

void ChemTable::update_number_of_analytical_species()
{
    auto count = std::count_if(model_.cbegin(), model_.cend(), [](LogKModel type){ return type == LogKModel::ANA; });
    size_analytical_ = static_cast<int>(count);
}

void ChemTable::update_logK_analytical(double T,double P)
{
    if(size_analytical_ == 0)   return;

    const double T_inv = 1.0 / T;
    const double T_inv2 = T_inv*T_inv;
    const double T2 = T * T;
    double loge = 0.43429448190325182765112891891661;

    const double P_ref = 1e5; // 1bar

    const int size_hkf = noRows_ - size_analytical_;
    for (int i = size_hkf; i < noRows_; ++i)
    {
        logK_[i] = a1_[i] + a2_[i] * T + a3_[i] * T_inv + a4_[i] * std::log10(T) + c1_[i] * T_inv2 + c2_[i]*T2;
        logK_[i] += (P-P_ref)*mol_volume_[i]/(PhysicalConstants::IdealGasConstant*T)*loge;
    }

}

void ChemTable::fix_hkf_units(ChemTable::Type type)
{
    constexpr double cal2J = UnitConversionFactors::cal2J_;
    [[maybe_unused]] constexpr double bar2Pa = UnitConversionFactors::bar2Pa_;
    constexpr double eta = UnitConversionFactors::eta_;

    if(type == ChemTable::Type::Basis)
    {
        const int size_hkf = noRows_ - size_analytical_;

        for (int i = 0; i < size_hkf; ++i)
        {
            deltaG_[i] *= cal2J;
            deltaH_[i] *= cal2J;
            S_[i] *= cal2J;
            // a1_[i] *= 1e-1*cal2J / bar2Pa; correction @ah 26/2 2025
            // a2_[i] *= 1e+2*cal2J; correction @ah 26/2 2025
            //a3_[i] *= cal2J / bar2Pa; correction @ah 26/2 2025
            //a4_[i] *= 1e+4*cal2J; correction @ah 26/2 2025
// ----- Free Energy Calc --------------
           // a1_[i] *= 1e-1 * cal2J;
           // a2_[i] *= 1e-4 * cal2J;

           // a3_[i] *= cal2J;
           // a4_[i] *= 1e-4 * cal2J;
// ------------------------------------- No bar to Pa
            a1_[i] *= 1e-1*cal2J ;
             a2_[i] *= 1e+2*cal2J;
            a3_[i] *= cal2J;
            a4_[i] *= 1e+4*cal2J;
// ----------------------------------------------
            c1_[i] *= cal2J;
            c2_[i] *= 1e+4*cal2J;
        }

        multiply_vector_by_value(mol_weight_, 1.0e-3);

        for (int i = 0; i < size_hkf; ++i)
        {
            if (omega_[i] == 0.0 && charge_[i] == 0.0)
            {
                re_[i] = 0.0;
            }
            else
            {
                re_[i] = charge_[i] * charge_[i] / (omega_[i] / eta + charge_[i] / 3.082);
            }
            omega_[i] *= 1.e+5*cal2J;
        }
    }
    else if(type == ChemTable::Type::Complex)
    {
        const int size_hkf = noRows_ - size_analytical_;

        for (int i = 0; i < size_hkf; ++i)
        {
            deltaG_[i] *= cal2J;
            deltaH_[i] *= cal2J;
            S_[i] *= cal2J;
 //           a1_[i] *= 1e-1*cal2J / bar2Pa; correction @ah 26/2 2025
 //           a2_[i] *= 1e+2*cal2J;
 //           a3_[i] *= cal2J / bar2Pa;
 //           a4_[i] *= 1e+4*cal2J; correction @ah 26/2 2025
// ----- Free Energy Calc --------------
//            a1_[i] *= 1e-1 * cal2J;
//            a2_[i] *= 1e-4 * cal2J;
//            a3_[i] *= cal2J;
//            a4_[i] *= 1e-4 * cal2J;
// -------------------------------------
            // ------------------------------------- No bar to Pa
            a1_[i] *= 1e-1 * cal2J;
            a2_[i] *= 1e+2 * cal2J;
            a3_[i] *= cal2J;
            a4_[i] *= 1e+4 * cal2J;
            // ----------------------------------------------
            c1_[i] *= cal2J;
            c2_[i] *= 1e+4*cal2J;
        }

        multiply_vector_by_value(mol_weight_, 1.0e-3);

        for (int i = 0; i < size_hkf; ++i)
        {
            if (omega_[i] == 0.0 && charge_[i] == 0.0)
            {
                re_[i] = 0.0;
            }
            else
            {
                re_[i] = charge_[i] * charge_[i] / (omega_[i] / eta + charge_[i] / 3.082);
            }
            omega_[i] *= 1.e+5*cal2J;
        }
    }
    else if(type == ChemTable::Type::Mineral)
    {
        const int size_hkf = noRows_ - size_analytical_;

        for (int i = 0; i < size_hkf; ++i)
        {
            deltaG_[i] *= cal2J;
            deltaH_[i] *= cal2J;
            S_[i] *= cal2J;
            a1_[i] *= cal2J;
            a2_[i] *= 1e-3*cal2J;
            a3_[i] *= 1e+5*cal2J;
        }
        multiply_vector_by_value(mol_volume_, 1.0e-6);
        multiply_vector_by_value(mol_weight_, 1.0e-3);

        // Finally, need to fix analytical (only for minerals)
        for (int k = size_hkf; k < noRows_; ++k)
        {   // logK_[i] = a1_[i] + a2_[i] * T + a3_[i] * T_inv + a4_[i] * logT + c1_[i] * T_inv2 + c2_[i]*T*T;
//            c1_[k]  = a7_[k];
//            a4_[k]  = a3_[k];
//            a3_[k]  = a2_[k];
//            a2_[k]  = a1_[k];
//            a1_[k]  = S_[k];
//            S_[k]   = 0.0;
            c2_[k]  = a7_[k];
            c1_[k]  = a3_[k];
            a4_[k]  = a2_[k];
            a3_[k]  = a1_[k];
            a2_[k]  = S_[k];
            a1_[k] = deltaH_[k];
            S_[k]   = 0.0;
            deltaH_[k] = 0.0;


        }
    }
}

void ChemTable::set_up_sparse_matrix()
{
    sparseM_.clear();
    sparseM_.resize(noRows_);

    for (int i = 0; i < noRows_; ++i)
    {
        for (int j = 0; j < noColumns_; ++j)
        {
            if (M_[i][j] != 0)
            {
                sparseM_[i].push_back({ j, M_[i][j] });
            }
        }
    }
}

void ChemTable::erase_row(std::size_t pos)
{
    // CHECK: ">" or ">=" ?
    if (row_name_.size() > pos)     remove_element_from_vector(row_name_, pos);
    if (nick_name_.size() > pos)    remove_element_from_vector(nick_name_, pos);
    if (abs_pos_.size() > pos)      remove_element_from_vector(abs_pos_, pos);
    if (log_m_.size() > pos)        remove_element_from_vector(log_m_, pos);
    if (log_a_.size() > pos)        remove_element_from_vector(log_a_, pos);
    if (log_g_.size() > pos)        remove_element_from_vector(log_g_, pos);
    if (log_QK_.size() > pos)       remove_element_from_vector(log_QK_, pos);
    if (delta_.size() > pos)        remove_element_from_vector(delta_, pos);
    if (logK_.size() > pos)         remove_element_from_vector(logK_, pos);
    if (dG_TP_.size() > pos)        remove_element_from_vector(dG_TP_, pos);
    if (a0_.size() > pos)           remove_element_from_vector(a0_, pos);
    if (a1_.size() > pos)           remove_element_from_vector(a1_, pos);
    if (a2_.size() > pos)           remove_element_from_vector(a2_, pos);
    if (a3_.size() > pos)           remove_element_from_vector(a3_, pos);
    if (a4_.size() > pos)           remove_element_from_vector(a4_, pos);
    if (a7_.size() > pos)           remove_element_from_vector(a7_, pos);
    if (c1_.size() > pos)           remove_element_from_vector(c1_, pos);
    if (c2_.size() > pos)           remove_element_from_vector(c2_, pos);
    if (deltaG_.size() > pos)       remove_element_from_vector(deltaG_, pos);
    if (deltaH_.size() > pos)       remove_element_from_vector(deltaH_, pos);
    if (S_.size() > pos)            remove_element_from_vector(S_, pos);
    if (omega_.size() > pos)        remove_element_from_vector(omega_, pos);
    if (re_.size() > pos)           remove_element_from_vector(re_, pos);
    if (charge_.size() > pos)       remove_element_from_vector(charge_, pos);
    if (scharge_.size() > pos)      remove_element_from_vector(scharge_, pos);
    if (mol_volume_.size() > pos)   remove_element_from_vector(mol_volume_, pos);
    if (mol_weight_.size() > pos)   remove_element_from_vector(mol_weight_, pos);
    if (log_af_.size() > pos)       remove_element_from_vector(log_af_, pos);
    if (type_.size() > pos)         remove_element_from_vector(type_, pos);
    if (model_.size() > pos)        remove_element_from_vector(model_, pos);
    if (T_.size() > pos)            remove_element_from_vector(T_, pos);
    if (M_.size() > pos)            remove_element_from_vector(M_, pos);
    --noRows_;
}

/** Erase table elements based on the names given as input. */
void ChemTable::erase_row_elements(const std::vector<std::string>& names)
{
    for(const auto& name: names)
    {
        int posi = get_row_index(name);
        if (posi > -1)
        {
            erase_row(posi);
        }
    }
}

void ChemTable::resize_col(int N, double default_val)
{
    col_name_.resize(N);

    for (int i = 0; i < noRows_; ++i)
    {
        M_[i].resize(N, default_val);
    }
    noColumns_ = N;
}

void ChemTable::resize_row(int N)
{
    // The first two are set properly outside, immediately after calling this function
    row_name_.resize(N, "");
    nick_name_.resize(N, "");

    if (!model_.empty())            model_.resize(N, LogKModel::HKF);
    if (!log_m_.empty())            log_m_.resize(N, 0.0);
    if (!log_a_.empty())            log_a_.resize(N, 0.0);
    if (!log_g_.empty())            log_g_.resize(N, 0.0);
    if (!log_af_.empty())           log_af_.resize(N, 0.0);
    if (!delta_.empty())            delta_.resize(N, 0.00);
    if (!T_.empty())                T_.resize(N, 0.0);
    if (!log_QK_.empty())           log_QK_.resize(N, 0.0);
    if (!mol_volume_.empty())       mol_volume_.resize(N, 0.0);
    if (!mol_weight_.empty())       mol_weight_.resize(N, 0.0);
    if (!a0_.empty())               a0_.resize(N, 0.0);
    if (!charge_.empty())           charge_.resize(N, 0.0);
    if (!scharge_.empty())          scharge_.resize(N, 0.0);
    if (!logK_.empty())             logK_.resize(N, 0.0);
    if (!abs_pos_.empty())          abs_pos_.resize(N, -1);
    if (!abs_pos_bas_.empty())      abs_pos_bas_.resize(N, -1);
    if (!type_.empty())             type_.resize(N, -1);
    if (!dG_TP_.empty())            dG_TP_.resize(N, 0.0);
    if (!deltaG_.empty())           deltaG_.resize(N, 0.0);
    if (!deltaH_.empty())           deltaH_.resize(N, 0.0);
    if (!S_.empty())                S_.resize(N, 0.0);
    if (!a1_.empty())               a1_.resize(N, 0.0);
    if (!a3_.empty())               a2_.resize(N, 0.0);
    if (!a3_.empty())               a3_.resize(N, 0.0);
    if (!a4_.empty())               a4_.resize(N, 0.0);
    if (!c1_.empty())               c1_.resize(N, 0.0);
    if (!c2_.empty())               c2_.resize(N, 0.0);
    if (!omega_.empty())            omega_.resize(N, 0.0);
    if (!re_.empty())               re_.resize(N, 0.0);

    M_.resize(N);
    for (int i = noRows_; i < N; ++i)
    {
        M_[i].resize(noColumns_, 0.0);
    }

    noRows_ = N;
}

/** Private helper function used when creating reduced ChemTables. */
void ChemTable::copy_data_to_new_table(ChemTable& CM) const
{
    copy_ChemTable_row(CM.T_, CM, T_);
    //
    copy_ChemTable_column(CM.model_, CM, model_);
    //
    copy_ChemTable_column(CM.delta_, CM, delta_);
    copy_ChemTable_column(CM.fugacity_, CM, fugacity_);
    copy_ChemTable_column(CM.log_QK_, CM, log_QK_);
    //
    copy_ChemTable_column(CM.type_ , CM, type_);
    copy_ChemTable_column(CM.logK_, CM, logK_);
    copy_ChemTable_column(CM.log_af_, CM, log_af_);
    copy_ChemTable_column(CM.a0_, CM, a0_);
    copy_ChemTable_column(CM.charge_, CM, charge_);
    copy_ChemTable_column(CM.scharge_, CM, scharge_);
    copy_ChemTable_column(CM.mol_volume_, CM, mol_volume_);
    copy_ChemTable_column(CM.mol_weight_, CM, mol_weight_);
    //
    copy_ChemTable_column(CM.deltaG_, CM, deltaG_);
    copy_ChemTable_column(CM.deltaH_, CM, deltaH_);
    copy_ChemTable_column(CM.S_, CM, S_);
    copy_ChemTable_column(CM.a1_, CM, a1_);
    copy_ChemTable_column(CM.a2_, CM, a2_);
    copy_ChemTable_column(CM.a3_, CM, a3_);
    copy_ChemTable_column(CM.a4_, CM, a4_);
    copy_ChemTable_column(CM.a7_, CM, a7_);
    copy_ChemTable_column(CM.c1_, CM, c1_);
    copy_ChemTable_column(CM.c2_, CM, c2_);
    copy_ChemTable_column(CM.omega_, CM, omega_);
    //
    copy_ChemTable_column(CM.re_, CM, re_);

    // gas EOS
    copy_ChemTable_column(CM.Tc_, CM, Tc_);
    copy_ChemTable_column(CM.Pc_, CM, Pc_);
    copy_ChemTable_column(CM.omega_acc_, CM, omega_acc_);
}

std::string ChemTable::get_reaction(const int i, int prec) const
{
    std::string reac;
    reac = row_name_[i] + " = ";
    for(int j=0;j<noColumns_;++j)
    {
        std::stringstream stream;
        if(prec>0)
            stream << std::fixed << std::setprecision(0) << M_[i][j];
        else
            stream << int(M_[i][j]);
        std::string ff = stream.str();
        if (M_[i][j] > 0)
        {
            reac += "+" + ff+ col_name_[j];
        }
        else if (M_[i][j] < 0)
        {
            reac += "-"+ ff + col_name_[j];
        }

        }
    return reac;
}
// sets new_value and returns old value
// returns tru if succesful

bool ChemTable::set_value(const std::string name_of_vector, const int index, const double new_value)
{
    if (!M_.empty() && !col_name_.empty())
    {
        if (col_name_.size() == M_[0].size()) // should set the value in matrix as well if the value is there as well
        {
            std::vector<std::string>::iterator it;
            for (it=col_name_.begin();it!=col_name_.end();++it)
            {
                if (*it == name_of_vector)
                {
                    int idx = it - col_name_.begin();
                    M_[index][idx] = new_value;
                }
            }
        }
    }
    auto pos = double_vector_map_.find(name_of_vector);
    if (pos != double_vector_map_.end()) //exists?
    {
        if (pos->second.size() > 0)
        {
           pos->second[index] = new_value;
            return true;
        }
        else
            return false;
    }
    else
        return false;
}

bool ChemTable::set_value(const std::string name_of_vector, const int index, const std::string new_value)
{
    auto pos = string_vector_map_.find(name_of_vector);
    if (pos != string_vector_map_.end()) //exists?
    {
        if (pos->second.size() > 0)
        {
            pos->second[index] = new_value;
            return true;
        }
        else
            return false;
    }
    else
        return false;
}

// find value in table and return true if found
bool ChemTable::get_value(const std::string name_of_vector, const int index, double &return_val)
{
    auto pos = double_vector_map_.find(name_of_vector);
    if (pos != double_vector_map_.end()) //exists?
    {
        if (pos->second.size()>0)
        {
            return_val = pos->second[index];
            return true;
        }
        else
            return false;
    }
    else
        return false;
}

// find value in table and return true if found
bool ChemTable::get_value(const std::string name_of_vector, const int index, std::string& return_val)
{
    auto pos = string_vector_map_.find(name_of_vector);
    if (pos != string_vector_map_.end()) //exists?
    {
        if (pos->second.size() > 0)
        {
            return_val = pos->second[index];
            return true;
        }
        else
            return false;
    }
    else
        return false;
}
