

#include "block_processing.h"

#include <list>

#include <utility>



void BlockProcess::SetCurBlock(uint64_t _cur_no_vec, uint8_t *cur_data)

{

    data = cur_data;
    cur_vec_len = params.vec_len;

    cur_no_vec = _cur_no_vec;

}

void BlockProcess::ProcessSquareBlock(uint32_t col_block_size, uint32_t col_vec_len, vector<uint32_t> &perm, vector<bool> &zeros, vector<bool> &copies, vector<uint32_t> &origin_of_copy, vector<uint8_t> &samples_indexes,bool permute)

{

    // Use per-column vec_len when tiling; fall back to global params for legacy path.
    cur_vec_len = col_vec_len ? col_vec_len : params.vec_len;
    // Clear per-block state to avoid cross-column contamination
    start = 0;

    if (permute)

    {

        permute_range_vec(0, cur_no_vec, perm,zeros, copies, origin_of_copy,samples_indexes);
    }

    else

    {

        perm.clear();

        perm.resize(col_block_size ? col_block_size : params.vec_len * 8, 0);

        for (int i = 0; i < (int)perm.size(); ++i)

            perm[i] = i;
    }
}


void BlockProcess::permute_range_vec(uint64_t id_start, uint64_t id_stop, vector<uint32_t> &v_perm,vector<bool> &zeros, vector<bool> &copies, vector<uint32_t> &origin_of_copy, vector<uint8_t> &samples_indexes)

{
   

    size_t n_h_samples = v_perm.size();

    // Use number of vectors in this block (rows) to size working buffers.
    uint64_t part_vec = id_stop - id_start;
    size_t max_no_vec_in_block = static_cast<size_t>(part_vec);
    
    comp_pos_copy_buf_.assign(max_no_vec_in_block, 0);

    // uint8_t* comp_samples_indexes = new uint8_t[max_no_vec_in_block*max_no_vec_in_block];
    sparse_matrix_cols_buf_.clear();
    if (sparse_matrix_cols_buf_.capacity() < max_no_vec_in_block * n_h_samples)
        sparse_matrix_cols_buf_.reserve(max_no_vec_in_block * n_h_samples);

    no_copy = 0;

    no_samples_index = 0;

    

    const uint32_t MC_ARRAY_SIZE = (max_no_vec_in_block+63) / 64;          

    struct mc_span_t
    {
        const uint64_t *p = nullptr;
        size_t n = 0;
        size_t size() const { return n; }
        const uint64_t *data() const { return p; }
        const uint64_t &operator[](size_t i) const { return p[i]; }
    };

    // One contiguous matrix: [n_h_samples][MC_ARRAY_SIZE]
    mc_vectors_buf_.assign(n_h_samples * MC_ARRAY_SIZE, 0);
    auto mc_at = [&](size_t i, size_t j) -> uint64_t& {
        return mc_vectors_buf_[i * MC_ARRAY_SIZE + j];
    };
    auto mc_row = [&](int i) -> mc_span_t {
        return mc_span_t{mc_vectors_buf_.data() + static_cast<size_t>(i) * MC_ARRAY_SIZE, MC_ARRAY_SIZE};
    };
    n_ones_buf_.assign(n_h_samples, 0);

    mc_ids_buf_.clear();

    for (int i = 0; i < (int)part_vec; ++i)

    {

        uint32_t id_cur = i + id_start;

        auto cur_vec = data + id_cur * cur_vec_len;

        bool empty = true;

        for (int j = 0; j < (int)cur_vec_len; ++j)

            if (cur_vec[j])

            {

                empty = false;

                break;
            }

        if (!empty)

            mc_ids_buf_.emplace_back(i);
        else{
            zeros[i] = 1;
            // cout<<"zi:"<<i<<endl;
        }
    }

    // random_shuffle(mc_ids.begin(), mc_ids.end());

    uint32_t mc_sort_size = mc_ids_buf_.size() >  max_no_vec_in_block ? max_no_vec_in_block  : mc_ids_buf_.size(); 

    // sort(mc_ids.begin(), mc_ids.begin() + mc_sort_size);
    // cout<<"mc_sort_size:"<<mc_sort_size<<endl;

    for (uint32_t i = 0; i < mc_sort_size; ++i)

    {

        uint32_t id_cur = mc_ids_buf_[i] + id_start;

        auto cur_vec = data + id_cur * cur_vec_len;

        uint32_t arr_id = i / 64;

        uint32_t arr_pos = i % 64;

        for (size_t j = 0; j < n_h_samples; ++j)

	        {

	            if (cur_vec[j / 8] & perm_lut8[j % 8])

	            {

	                mc_at(j, arr_id) |= perm_lut64[arr_pos];

	                ++n_ones_buf_[j];
	            }
	        }
	    }
    // Keep capacity for reuse across blocks (avoids allocator churn).
    // Determine the best permutation

    perm_buf_.clear();
    perm_buf_.reserve(n_h_samples + 2);

    uint64_t cost_perm = 0;

    

    list<pair<int, int>> density_list;

    density_list.emplace_back(make_pair(0, n_ones_buf_[0]));
    perm_buf_.emplace_back(0);
    auto p = density_list.begin();

    for (int i = 1; i < (int)n_ones_buf_.size(); ++i)

        density_list.emplace_back(make_pair(i, n_ones_buf_[i]));

    // Insert two guards into list

    int huge_val = 1 << 28;

    density_list.emplace_back(make_pair(-1, -2 * huge_val));

    density_list.emplace_back(make_pair(-1, 2 * huge_val));

    density_list.sort([](pair<int, int> &x, pair<int, int> &y)
                      { return x.second < y.second; });
        
    // auto p = density_list.begin();
    // ++p;
    // perm.emplace_back(p->first);
    // cout<<p->first<<":"<<p->second<<endl;
    size_t n_cur_samples = n_h_samples;

    new_mc_vectors_buf_.resize(n_h_samples * MC_ARRAY_SIZE);
    auto new_mc_at = [&](size_t i, size_t j) -> uint64_t& {
        return new_mc_vectors_buf_[i * MC_ARRAY_SIZE + j];
    };

    new_xor_buf_.resize(MC_ARRAY_SIZE);
    best_new_xor_buf_.resize(MC_ARRAY_SIZE);
    
    while (n_cur_samples)

    {
        
	        for(size_t t = 0 ; t < MC_ARRAY_SIZE; t++)        
	        {
	            if((n_h_samples-n_cur_samples) %8 == 0)
	            
	                best_new_xor_buf_[t] = mc_at(static_cast<size_t>(p->first), t);
	            
	            new_mc_at(n_h_samples-n_cur_samples, t)=best_new_xor_buf_[t];

	        }     

        uint64_t best_cost = huge_val;

        auto best_p = p;


        // Look for the most similar vector to *p

        // Starts from p and moves up and down on the list ordered according to the number of ones

        // This limits the number of vector pairs that must be evaluated

        auto p_down = p;

        auto p_up = p;

        --p_down;

        ++p_up;

        uint64_t dif_up = abs(p->second - p_up->second);

        uint64_t dif_down = abs(p->second - p_down->second);

        while (true)
        {
            uint64_t min_dif = min(dif_up, dif_down);

            if (min_dif >= best_cost)

                break;

	            if (dif_up < dif_down)

	            {

		                uint64_t cost = bit_cost(mc_row(p->first), mc_row(p_up->first), best_cost, new_xor_buf_.data());
	                if (cost < best_cost)

	                {

	                    best_cost = cost;

                    best_p = p_up;
                   
	                    memcpy(best_new_xor_buf_.data(),new_xor_buf_.data(),MC_ARRAY_SIZE*sizeof(uint64_t));

                    if (best_cost == 0)

                        break;
                }

                ++p_up;

                dif_up = abs(p->second - p_up->second);
            }

	            else
	            {
		                uint64_t cost = bit_cost(mc_row(p->first), mc_row(p_down->first), best_cost, new_xor_buf_.data());
	                
	                if (cost < best_cost)

	                {

                    best_cost = cost;

                    best_p = p_down;

	                    memcpy(best_new_xor_buf_.data(),new_xor_buf_.data(),MC_ARRAY_SIZE*sizeof(uint64_t));

                    if (best_cost == 0)

                        break;
                }

                --p_down;

                dif_down = abs(p->second - p_down->second);
            }
        }

        perm_buf_.emplace_back(best_p->first);

        density_list.erase(p);

        p = best_p;

        cost_perm += best_cost;

        --n_cur_samples;
    }

    v_perm.clear();

    v_perm.resize(n_h_samples, 0);

    for (size_t i = 0; i < n_h_samples; ++i)
    {
        // v_perm[perm[i]] = i;
        v_perm[i] = perm_buf_[i];

    }
    old_vec_storage_.resize(cur_vec_len);
    uint8_t *old_vec = old_vec_storage_.data();

    size_t mc_tr=0;

    const uint32_t NO_COPY = std::numeric_limits<uint32_t>::max();
    copy_from_buf_.assign(static_cast<size_t>(part_vec), NO_COPY);

    for (int i = 0; i < (int)part_vec; ++i)

    {
        

        bool empty = true;

        auto cur_arr_id = (i-mc_tr)/64;

        auto cur_arr_pos = (i-mc_tr)%64; 

        uint32_t id_cur = i + id_start;

        auto new_vec = data + id_cur * cur_vec_len;

        memcpy(old_vec, new_vec, cur_vec_len);

        fill_n(new_vec, cur_vec_len, 0);

        for (int j = 0; j < (int)cur_vec_len; ++j){

            if (old_vec[j])

            {
                empty = false;

                break;
            }
        }

        if(empty){
            mc_tr++;
            continue;
        }
        else {
            if(!copies[i]){
                for(size_t next_copy = i+1 ; next_copy < part_vec && next_copy < i+params.max_replication_depth;next_copy++){

                    auto next_vec = data+ next_copy * cur_vec_len;

                    if(memcmp(old_vec,next_vec,cur_vec_len) == 0)
                    {
                        if (!copies[next_copy])
                        {
                            copies[next_copy] = 1;
                            if (copy_from_buf_[next_copy] == NO_COPY)
                                copy_from_buf_[next_copy] = static_cast<uint32_t>(i);
                        }
                    // comp_pos_copy[no_copy++] = i;
                    // cout<<"i:"<<i<<endl;
   
                    }
                }
                uint32_t prev_index = 0;
                for(size_t t = 1 ; t<= n_h_samples;t++){
                
	                    if(new_mc_at(t-1,cur_arr_id)&perm_lut64[cur_arr_pos]){
                        
                        sparse_matrix_cols_buf_.emplace_back(static_cast<uint32_t>(t) - prev_index);
                        prev_index = static_cast<uint32_t>(t) ;

                    }
                }
                sparse_matrix_cols_buf_.emplace_back(0);
  
                // uint8_t first_prev_index = 0;
                // uint8_t second_prev_index =0;
                // for(size_t t = 0 ;t<n_h_samples;t++){
                
                //     if(new_mc_vectors[t][cur_arr_id]&perm_lut64[cur_arr_pos]){
                //     // cout<<t<<" ";
                //         comp_samples_indexes[no_samples_index++] = (t & 0xffu)-first_prev_index;
                        
                //         first_prev_index = (t & 0xffu) ;

                //         comp_samples_indexes[no_samples_index++] = ((t >>8) & 0xffu)-second_prev_index;

                //         second_prev_index = ((t >>8) & 0xffu) ;

                //     }
                // }
            

                // comp_samples_indexes[no_samples_index++] = 0;
                // comp_samples_indexes[no_samples_index++] = 0;

            }
            else

                continue;
            
                

        }
    }
    for (size_t i = 0; i < copy_from_buf_.size(); ++i)
    {
        if (copy_from_buf_[i] != NO_COPY)
            comp_pos_copy_buf_[no_copy++] = copy_from_buf_[i];
    }
    
    origin_of_copy.assign(comp_pos_copy_buf_.begin(), comp_pos_copy_buf_.begin() + no_copy);
   
    vint_code::EncodeArray(sparse_matrix_cols_buf_, samples_indexes);
    sparse_matrix_cols_buf_.clear();

    // samples_indexes.resize(no_samples_index);

    // memcpy(&samples_indexes[0], comp_samples_indexes, no_samples_index* sizeof(uint8_t));

    // block_index_size = no_samples_index;

    //_copy_num = no_copy;

    // new_mc_vectors and comp_pos_copy are std::vector-backed; no manual cleanup needed.
    // delete [] comp_samples_indexes;

}
void BlockProcess::ProcessLastBlock(vector<bool> &zeros, vector<bool> &copies, vector<uint32_t> &origin_of_copy, vector<uint8_t> &samples_indexes )
{
    const uint32_t NO_COPY = std::numeric_limits<uint32_t>::max();
    copy_from_buf_.assign(static_cast<size_t>(cur_no_vec), NO_COPY);
    comp_pos_copy_buf_.assign(static_cast<size_t>(cur_no_vec), 0);
    // encode_byte_num = (uint32_t)log2(params.n_samples*2-1)/8+1;
    // uint8_t* comp_samples_indexes = new uint8_t[cur_no_vec*cur_vec_len*8*2];
    sparse_matrix_cols_buf_.clear();
    if (sparse_matrix_cols_buf_.capacity() < cur_no_vec * cur_vec_len * 8 * 2)
        sparse_matrix_cols_buf_.reserve(cur_no_vec * cur_vec_len * 8 * 2);
    no_copy = 0;
    no_samples_index = 0;

    for (int i = 0; i < (int)cur_no_vec; ++i)

    {
        
        bool empty = true;

        auto cur_vec = data + i * cur_vec_len;

        for (int j = 0; j < (int)cur_vec_len; ++j){

            if (cur_vec[j])

            {
                empty = false;

                break;
            }
        }
        if(empty){
            zeros[i] =1;
            // cout<<"zi1:"<<i<<endl;
            continue;
        }
        else {
            if(!copies[i]){
                for(size_t next_copy = i+1 ; next_copy < cur_no_vec && next_copy < i+params.max_replication_depth;next_copy++){

                    auto next_vec = data+ next_copy * cur_vec_len;

                    if(memcmp(cur_vec,next_vec,cur_vec_len) == 0)
                    {
                        if (!copies[next_copy])
                        {
                            copies[next_copy] = 1;
                            if (copy_from_buf_[next_copy] == NO_COPY)
                                copy_from_buf_[next_copy] = static_cast<uint32_t>(i);
                        }
                    }
                }
                uint32_t prev_index = 0;
                for(size_t t = 1 ; t<= cur_vec_len*8;t++){
                    auto cur_arr_pos = (t-1)%8;
                    if(cur_vec[(t-1)/8]&perm_lut8[cur_arr_pos]){
                         
                        sparse_matrix_cols_buf_.emplace_back(static_cast<uint32_t>(t) - prev_index);
                        prev_index = static_cast<uint32_t>(t) ;

                    }
                }
                sparse_matrix_cols_buf_.emplace_back(0);

            //     uint8_t first_prev_index = 0;            
            //     uint8_t second_prev_index =0;            

            //     for(size_t t = 0 ;t<params.vec_len*8;t++){
            //     // cout<<new_mc_vectors[t][cur_arr_id]<<endl;
                    
            //         auto cur_arr_pos = t%8; 

            //         if(cur_vec[t/8] & perm_lut8[cur_arr_pos]){

            //             // // comp_samples_indexes[no_samples_index++] = t & 0xffu;  
                        
            //             // // comp_samples_indexes[no_samples_index++] = (t >>8) & 0xffu;

            //             comp_samples_indexes[no_samples_index++] = (t & 0xffu)-first_prev_index;  
                        
            //             first_prev_index = (t & 0xffu) ;

            //             comp_samples_indexes[no_samples_index++] = ((t >>8) & 0xffu)-second_prev_index;

            //             second_prev_index = ((t >>8) & 0xffu) ;

                        
            //         }
            // //    cout<<size_t(new_vec[t/8])<<" ";    
            //     }
            //     comp_samples_indexes[no_samples_index++] = 0;
            //     comp_samples_indexes[no_samples_index++] = 0;
                
            }
            else

                continue;

                

        }
    }

    for (size_t i = 0; i < copy_from_buf_.size(); ++i)
    {
        if (copy_from_buf_[i] != NO_COPY)
            comp_pos_copy_buf_[no_copy++] = copy_from_buf_[i];
    }

    origin_of_copy.assign(comp_pos_copy_buf_.begin(), comp_pos_copy_buf_.begin() + no_copy);

    vint_code::EncodeArray(sparse_matrix_cols_buf_, samples_indexes);
    sparse_matrix_cols_buf_.clear();
    // samples_indexes.resize(no_samples_index); 

    // memcpy(&samples_indexes[0], comp_samples_indexes, no_samples_index* sizeof(uint8_t));

    // block_index_size = no_samples_index;

   // _copy_num = no_copy;
    //  cout<<"samples_indexs:";
    // for(size_t i = 0;i<no_samples_index;i++)
        // cout<<block_index_size<<endl;
    // delete [] comp_samples_indexes;
   

}
void BlockProcess::ProcessVariant(vector<uint32_t> &perm, vector<variant_desc_t> &v_vcf_data_io)
{

    if (v_vcf_data_io.size() == perm.size())
    {
        get_perm(perm, perm.size(), v_vcf_data_io);
    }
}
inline void BlockProcess::get_perm(const vector<uint32_t> &perm, int n, vector<variant_desc_t> &v_vcf_data_compress)
{
    if (n <= 0)
        return;
    if (perm_pos_buf_.size() < static_cast<size_t>(n))
        perm_pos_buf_.resize(n);
    if (perm_ref_buf_.size() < static_cast<size_t>(n))
        perm_ref_buf_.resize(n);

    for (int j = 0; j < n; j++)
    {
        int i = perm[j];
        perm_pos_buf_[i] = v_vcf_data_compress[j].pos;
        perm_ref_buf_[i] = std::move(v_vcf_data_compress[j].ref);
    }

    for (int j = 0; j < n; j++)
    {
        v_vcf_data_compress[j].pos = perm_pos_buf_[j];
        v_vcf_data_compress[j].ref = std::move(perm_ref_buf_[j]);
    }
}
void BlockProcess::AddGtBlock(fixed_field_block &_fixed_field_block_io, vector<bool> &_all_zeros, vector<bool> &_all_copies, vector<uint32_t> &_comp_pos_copy,
                              vector<bool> &_zeros_only, vector<bool> &_copies, vector<uint32_t> &_origin_of_copy, vector<uint8_t> &_samples_indexes)
{

    for (size_t i = 0; i < _origin_of_copy.size(); i++)
    {
        _origin_of_copy[i] += _all_zeros.size();
    }

    _all_zeros.insert(_all_zeros.end(), _zeros_only.begin(), _zeros_only.end());
    _all_copies.insert(_all_copies.end(), _copies.begin(), _copies.end());
    _comp_pos_copy.insert(_comp_pos_copy.end(), _origin_of_copy.begin(), _origin_of_copy.end());
    if (!_samples_indexes.empty())
        _fixed_field_block_io.gt_block.reserve(_fixed_field_block_io.gt_block.size() + _samples_indexes.size());
    _fixed_field_block_io.gt_block.insert(_fixed_field_block_io.gt_block.end(), _samples_indexes.begin(), _samples_indexes.end());
}

void BlockProcess::AddGtIndexBlock(std::vector<uint8_t> &gt_block, vector<bool> &_all_zeros, vector<bool> &_all_copies, vector<uint32_t> &_comp_pos_copy,
                                  vector<bool> &_zeros_only, vector<bool> &_copies, vector<uint32_t> &_origin_of_copy, vector<uint8_t> &_samples_indexes)
{
    for (size_t i = 0; i < _origin_of_copy.size(); i++)
    {
        _origin_of_copy[i] += _all_zeros.size();
    }

    _all_zeros.insert(_all_zeros.end(), _zeros_only.begin(), _zeros_only.end());
    _all_copies.insert(_all_copies.end(), _copies.begin(), _copies.end());
    _comp_pos_copy.insert(_comp_pos_copy.end(), _origin_of_copy.begin(), _origin_of_copy.end());
    if (!_samples_indexes.empty())
        gt_block.reserve(gt_block.size() + _samples_indexes.size());
    gt_block.insert(gt_block.end(), _samples_indexes.begin(), _samples_indexes.end());
}

void BlockProcess::addFixedFieldsBlock(fixed_field_block &_fixed_field_block_io, vector<variant_desc_t> &_v_vcf_data_io, int64_t &prev_pos)
{
    _fixed_field_block_io.no_variants += static_cast<uint32_t>(_v_vcf_data_io.size());
    for (size_t i = 0; i < _v_vcf_data_io.size(); ++i)
    {
        const variant_desc_t &desc = _v_vcf_data_io[i];
        append_str(_fixed_field_block_io.chrom, desc.chrom);
        append_str(_fixed_field_block_io.id, desc.id);
        append_str(_fixed_field_block_io.alt, desc.alt);
        append_str(_fixed_field_block_io.qual, desc.qual);
        append(_fixed_field_block_io.pos, (int64_t)(desc.pos - prev_pos));
        prev_pos = desc.pos;
        append_str(_fixed_field_block_io.ref, desc.ref);
    }
}
void BlockProcess::addSortFieldBlock(fixed_field_block &_fixed_field_block_io,vector<bool> &_all_zeros,vector<bool> &_all_copies,vector<uint32_t> &_comp_pos_copy,
    
    vector<bool> &_zeros_only, vector<bool> &_copies, vector<uint32_t> &_origin_of_copy,vector<uint8_t> &_samples_indexes,vector<variant_desc_t> &_v_vcf_data_io,int64_t &prev_pos)
{         

    AddGtBlock(_fixed_field_block_io, _all_zeros, _all_copies, _comp_pos_copy, _zeros_only, _copies, _origin_of_copy, _samples_indexes);
    addFixedFieldsBlock(_fixed_field_block_io, _v_vcf_data_io, prev_pos);
}
// void BlockProcess::addSortFieldBlock(sort_field_block &sort_fixed_field_block_io,vector<bool> &_all_zeros,vector<bool> &_all_copies,vector<uint32_t> &_comp_pos_copy,
//     vector<bool> &_zeros_only, vector<bool> &_copies, vector<uint32_t> &_origin_of_copy,vector<uint8_t> &_samples_indexes,FieldsPackage &fields_pck,int64_t &prev_pos){
//     unique_lock<mutex> lck(mtx_v_part1);
// 	cv_v_part1.wait(lck, [&, this] {
// 		return cur_block_id == fields_pck.block_id;
// 		});    
//     _all_zeros.insert(_all_zeros.end(),_zeros_only.begin(),_zeros_only.end());
//     _all_copies.insert(_all_copies.end(),_copies.begin(),_copies.end());
//     for(size_t i = 0 ;i < _origin_of_copy.size(); i++)
//     {
//         _origin_of_copy[i] += start;
//     }
//     _comp_pos_copy.insert(_comp_pos_copy.end(),_origin_of_copy.begin(),_origin_of_copy.end());

//     sort_fixed_field_block_io.gt_block.insert(sort_fixed_field_block_io.gt_block.end(), _samples_indexes.begin(), _samples_indexes.end());
//     start += _zeros_only.size();
//     variant_desc_t desc;
//     for (size_t i = 0; i < fields_pck.v_vcf_data_compress.size(); ++i)
//     {
//         desc = fields_pck.v_vcf_data_compress[i];
//         append(sort_fixed_field_block_io.pos, (int64_t)(desc.pos - prev_pos));
// 	    prev_pos = desc.pos;
// 	    append_str(sort_fixed_field_block_io.ref, desc.ref);
//     }
//     lock_guard<mutex> lck(mtx_v_part1);
// 	++cur_block_id;
// 	cv_v_part1.notify_all();
    
// }
// void BlockProcess::ProcessLastPerm(vector<uint32_t> &perm,vector<vector<uint8_t>> &_vint_last_perm){
//     vector<uint8_t> vint_perm;
//     vint_perm = vint_code::EncodeArray(perm);
//     _vint_last_perm.emplace_back(vint_perm);
// }
