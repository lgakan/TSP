#include "TSP.hpp"

#include <algorithm>
#include <stack>
#include <optional>
#include <iostream>

std::ostream& operator<<(std::ostream& os, const CostMatrix& cm) {
    for (std::size_t r = 0; r < cm.size(); ++r) {
        for (std::size_t c = 0; c < cm.size(); ++c) {
            const auto& elem = cm[r][c];
            os << (is_inf(elem) ? "INF" : std::to_string(elem)) << " ";
        }
        os << "\n";
    }
    os << std::endl;

    return os;
}

/* PART 1 */

/**
 * Create path from unsorted path and last 2x2 cost matrix.
 * @return The vector of consecutive vertex.
 */
// TODO: Implement it!
path_t StageState::get_path() {

    // Adding last 2 vertex to list
    for(std::size_t row = 0; row < matrix_.size(); row++){
        for(std::size_t col = 0; col < matrix_.size(); col++){
            if(matrix_[row][col] != INF){
                append_to_path(vertex_t(row, col));
            }
        }
    }
    // Creating sorted path
    path_t sorted_path;
    //Adding first vertex from unsorted path as first sorted vertex
    sorted_path.push_back(unsorted_path_[0].col);

    //Adding new sorted vertexes
    while(sorted_path.size() < matrix_.size()){
        for(auto i : unsorted_path_){
            if(sorted_path.back() == i.row){
                sorted_path.push_back(i.col);
            }
        }
    }

    //Starting indexing from 1
    for(auto& i : sorted_path){
        i = i+1;
    }

    return sorted_path;
}

/**
 * Get minimum values from each row and returns them.
 * @return Vector of minimum values in row.
 */
std::vector<cost_t> CostMatrix::get_min_values_in_rows() const{
    std::vector<cost_t> min_values = {};

     for(std::size_t row = 0; row < matrix_.size(); row++){
         cost_t min_single_row_value = matrix_[row][0];
         for(std::size_t col = 0; col < matrix_[row].size(); col++){
             if(matrix_[row][col] < min_single_row_value){
                 min_single_row_value = matrix_[row][col];
             }
         }
         min_values.push_back(min_single_row_value);
     }
    return min_values;
}

/**
 * Reduce rows so that in each row at least one zero value is present.
 * @return Sum of values reduced in rows.
 */
cost_t CostMatrix::reduce_rows() {
    cost_t sum_of_min_values = 0;
    std::vector<cost_t> min_values = CostMatrix(matrix_).get_min_values_in_rows();
    for(std::size_t row = 0; row < matrix_.size(); row++){
        for(std::size_t col = 0; col < matrix_[row].size(); col++){
            if(matrix_[row][col] != INF){
                matrix_[row][col] -= min_values[row];
            }

        }
        if(min_values[row] != INF){
            sum_of_min_values += min_values[row];
        }
    }
    return sum_of_min_values;
}

/**
 * Get minimum values from each column and returns them.
 * @return Vector of minimum values in columns.
 */
std::vector<cost_t> CostMatrix::get_min_values_in_cols() const {
    std::vector<cost_t> min_values = {};

    for(std::size_t col = 0; col < matrix_.size(); col++){
        cost_t min_single_row_value = matrix_[0][col];
        for(std::size_t row = 0; row < matrix_[col].size(); row++){
            if(matrix_[row][col] < min_single_row_value){
                min_single_row_value = matrix_[row][col];
            }
        }
        min_values.push_back(min_single_row_value);
    }
    return min_values;
}

/**
 * Reduces rows so that in each column at least one zero value is present.
 * @return Sum of values reduced in columns.
 */
cost_t CostMatrix::reduce_cols() {
    cost_t sum_of_min_values = 0;
    std::vector<cost_t> min_values = CostMatrix(matrix_).get_min_values_in_cols();

    for(std::size_t col = 0; col < matrix_.size(); col++){
        for(std::size_t row = 0; row < matrix_[col].size(); row++){
            if(matrix_[row][col] != INF){
                matrix_[row][col] -= min_values[col];
            }
        }
        if(min_values[col] != INF){
            sum_of_min_values += min_values[col];
        }
    }
    return sum_of_min_values;
}

/**
 * Get the cost of not visiting the vertex_t (@see: get_new_vertex())
 * @param row
 * @param col
 * @return The sum of minimal values in row and col, excluding the intersection value.
 */
cost_t CostMatrix::get_vertex_cost(std::size_t row, std::size_t col) const {
    cost_t row_min_value = INF;
    cost_t col_min_value = INF;

    for(std::size_t c = 0; c < matrix_.size(); c++){
        if(c == col){
            continue;
        }
        if(matrix_[row][c] < row_min_value){
            row_min_value = matrix_[row][c];
        }
    }

    for(std::size_t r = 0; r < matrix_.size(); r++){
        if(r== row){
            continue;
        }
        if(matrix_[r][col] < col_min_value){
            col_min_value = matrix_[r][col];
        }
    }
    return col_min_value + row_min_value;
}

/* PART 2 */

/**
 * Choose next vertex to visit:
 * - Look for vertex_t (pair row and column) with value 0 in the current cost matrix.
 * - Get the vertex_t cost (calls get_vertex_cost()).
 * - Choose the vertex_t with maximum cost and returns it.
 * @return The coordinates of the next vertex.
 */
NewVertex StageState::choose_new_vertex() {
    NewVertex max_cost_vertex = NewVertex(vertex_t(0, 0), -1);

    for(std::size_t row = 0; row < matrix_.size(); row++){
        for(std::size_t col = 0; col < matrix_[row].size(); col++){
            if(matrix_[row][col] == 0){
                vertex_t vertex = vertex_t(row, col);
                cost_t cost = matrix_.get_vertex_cost(vertex.row, vertex.col);
                if(cost > max_cost_vertex.cost){
                    max_cost_vertex.coordinates = vertex;
                    max_cost_vertex.cost = cost;
                }
            }
        }
    }
    return max_cost_vertex;
}

/**
 * Update the cost matrix with the new vertex.
 * @param new_vertex
 */
// TODO: Implement it! CYKLE
void StageState::update_cost_matrix(vertex_t new_vertex) {

    //Size of matrix reduction
    for(std::size_t row = 0; row < matrix_.size(); row++){
        for(std::size_t col = 0; col < matrix_[row].size(); col++){
            if(row == new_vertex.row or col == new_vertex.col){
                matrix_[row][col] = INF;
            }
        }
    }

    // if we get 1,3 we can't go 3,1
    matrix_[new_vertex.col][new_vertex.row] = INF;

    // Cycles
    std::vector<std::size_t> rows;
    std::vector<std::size_t> cols;

    for(auto i : unsorted_path_){
        rows.push_back(i.row);
        cols.push_back(i.col);
    }
    bool is_unique = true;
    std::vector<std::size_t> checking;
    size_t row_to_delete;
    size_t col_to_delete;

    // taking unique val
    for(std::size_t x = 0; x < rows.size(); x++){
        is_unique = true;

        for(std::size_t y = 0; y < cols.size(); y++){
            if(rows[x] == cols[y]){
                is_unique = false;
            }
        }

        if(is_unique){
            checking.push_back(rows[x]);
            row_to_delete = rows[x];
        }

        is_unique = true;
        for(std::size_t y = 0; y < cols.size(); y++){
            if(cols[x] == rows[y]){
                is_unique = false;
            }
        }
        if(is_unique){
            col_to_delete = cols[x];
        }
    }
    if(checking.size() == 1){
        matrix_[col_to_delete][row_to_delete] = INF;
    }

}

/**
 * Reduce the cost matrix.
 * @return The sum of reduced values.
 */
cost_t StageState::reduce_cost_matrix() {
    return matrix_.reduce_rows() + matrix_.reduce_cols();
}

/**
 * Given the optimal path, return the optimal cost.
 * @param optimal_path
 * @param m
 * @return Cost of the path.
 */
cost_t get_optimal_cost(const path_t& optimal_path, const cost_matrix_t& m) {
    cost_t cost = 0;

    for (std::size_t idx = 1; idx < optimal_path.size(); ++idx) {
        cost += m[optimal_path[idx - 1] - 1][optimal_path[idx] - 1];
    }

    // Add the cost of returning from the last city to the initial one.
    cost += m[optimal_path[optimal_path.size() - 1] - 1][optimal_path[0] - 1];

    return cost;
}


/**
 * Create the right branch matrix with the chosen vertex forbidden and the new lower bound.
 * @param m
 * @param v
 * @param lb
 * @return New branch.
 */
StageState create_right_branch_matrix(cost_matrix_t m, vertex_t v, cost_t lb) {
    CostMatrix cm(m);
    cm[v.row][v.col] = INF;
    return StageState(cm, {}, lb);
}

/**
 * Retain only optimal ones (from all possible ones).
 * @param solutions
 * @return Vector of optimal solutions.
 */
tsp_solutions_t filter_solutions(tsp_solutions_t solutions) {
    cost_t optimal_cost = INF;
    for (const auto& s : solutions) {
        optimal_cost = (s.lower_bound < optimal_cost) ? s.lower_bound : optimal_cost;
    }

    tsp_solutions_t optimal_solutions;
    std::copy_if(solutions.begin(), solutions.end(),
                 std::back_inserter(optimal_solutions),
                 [&optimal_cost](const tsp_solution_t& s) { return s.lower_bound == optimal_cost; }
    );

    return optimal_solutions;
}

/**
 * Solve the TSP.
 * @param cm The cost matrix.
 * @return A list of optimal solutions.
 */
tsp_solutions_t solve_tsp(const cost_matrix_t& cm) {

    StageState left_branch(cm);

    // The branch & bound tree.
    std::stack<StageState> tree_lifo;

    // The number of levels determines the number of steps before obtaining
    // a 2x2 matrix.
    std::size_t n_levels = cm.size() - 2;

    tree_lifo.push(left_branch);   // Use the first cost matrix as the root.

    cost_t best_lb = INF;
    tsp_solutions_t solutions;

    while (!tree_lifo.empty()) {

        left_branch = tree_lifo.top();
        tree_lifo.pop();

        while (left_branch.get_level() != n_levels && left_branch.get_lower_bound() <= best_lb) {
            // Repeat until a 2x2 matrix is obtained or the lower bound is too high...

            if (left_branch.get_level() == 0) {
                left_branch.reset_lower_bound();
            }

            // 1. Reduce the matrix in rows and columns.
//            cost_t new_cost = 0; // @TODO (KROK 1)
            cost_t new_cost = left_branch.reduce_cost_matrix();

            // 2. Update the lower bound and check the break condition.
            left_branch.update_lower_bound(new_cost);
            if (left_branch.get_lower_bound() > best_lb) {
                break;
            }

            // 3. Get new vertex and the cost of not choosing it.
            // @TODO (KROK 2)
            NewVertex new_vertex = left_branch.choose_new_vertex();

            // 4. @TODO Update the path - use append_to_path method.
            left_branch.append_to_path(new_vertex.coordinates);
            // 5. @TODO (KROK 3) Update the cost matrix of the left branch.
            left_branch.update_cost_matrix(new_vertex.coordinates);

            // 6. Update the right branch and push it to the LIFO.
            cost_t new_lower_bound = left_branch.get_lower_bound() + new_vertex.cost;
            tree_lifo.push(create_right_branch_matrix(cm, new_vertex.coordinates,
                                                      new_lower_bound));
        }

        if (left_branch.get_lower_bound() <= best_lb) {
            // If the new solution is at least as good as the previous one,
            // save its lower bound and its path.
            best_lb = left_branch.get_lower_bound();
            path_t new_path = left_branch.get_path();
            solutions.push_back({get_optimal_cost(new_path, cm), new_path});
        }
    }

    return filter_solutions(solutions); // Filter solutions to find only optimal ones.
}


