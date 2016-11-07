#theory dl {
    term{};
    constant {
    - : 1, unary};
    diff_term {
    - : 1, binary, left};
    &diff/0 : diff_term, {<=}, constant, any;
    &show_assignment/0 : term, directive
}.
