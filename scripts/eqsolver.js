// Equation Solver Module
const EqSolver = {
    equations: [],
    flags: {
        complex: false,
        integer: false,
        positive: false
    },

    reset() {
        this.equations = [];
        this.flags = {
            complex: false,
            integer: false,
            positive: false
        };
    },

    parseEquation(eqStr) {
        // Parse equation like "3x + y = 8" or "3x - y = 8 - y"
        const sides = eqStr.split('=');
        if (sides.length !== 2) {
            throw new Error('Invalid equation format. Use format: 3x + y = 8');
        }

        const lhs = sides[0].trim();
        const rhs = sides[1].trim();

        // Parse both sides as expressions with terms
        const lhsTerms = this.parseExpression(lhs);
        const rhsTerms = this.parseExpression(rhs);

        // Move all terms to left side (lhs - rhs = 0)
        const terms = { ...lhsTerms };
        
        for (const [variable, coef] of Object.entries(rhsTerms)) {
            if (variable === 'constant') {
                terms.constant = (terms.constant || 0) - coef;
            } else {
                terms[variable] = (terms[variable] || 0) - coef;
            }
        }

        // Extract constant term
        const constant = -(terms.constant || 0);
        delete terms.constant;

        if (Object.keys(terms).length === 0) {
            throw new Error('No variables found in equation');
        }

        return { terms, constant };
    },

    parseExpression(expr) {
        // Parse an expression and return terms with their coefficients
        const terms = { constant: 0 };
        
        // Remove all spaces first
        expr = expr.replace(/\s+/g, '');
        
        // Split into tokens: numbers, variables, operators
        // Match patterns like: +3x, -2y, +5, -10, x, y
        const pattern = /([+-])?(\d*\.?\d+)?([a-z])?/gi;
        const matches = expr.matchAll(pattern);
        
        let hasContent = false;
        
        for (const match of matches) {
            const fullMatch = match[0];
            if (!fullMatch) continue;
            
            const sign = match[1] || '+';
            const number = match[2];
            const variable = match[3] ? match[3].toLowerCase() : null;
            
            if (variable) {
                // Variable term like "3x" or "x"
                let coef = number ? parseFloat(number) : 1;
                if (sign === '-') coef = -coef;
                
                terms[variable] = (terms[variable] || 0) + coef;
                hasContent = true;
            } else if (number) {
                // Constant term like "5" or "-10"
                let value = parseFloat(number);
                if (sign === '-') value = -value;
                
                terms.constant += value;
                hasContent = true;
            }
        }

        if (!hasContent) {
            throw new Error('Invalid expression format');
        }

        return terms;
    },

    addEquation(eqStr) {
        const parsed = this.parseEquation(eqStr);
        this.equations.push(parsed);
        return this.equations.length - 1;
    },

    showEquations() {
        return this.equations.map((eq, idx) => {
            const terms = Object.entries(eq.terms)
                .map(([v, c]) => {
                    if (c === 0) return '';
                    const coefStr = c === 1 ? '' : (c === -1 ? '-' : c);
                    return `${c > 0 && idx > 0 ? '+' : ''}${coefStr}${v}`;
                })
                .filter(t => t)
                .join(' ');
            return `[${idx}] ${terms} = ${eq.constant}`;
        });
    },

    deleteEquations(indices) {
        // Sort indices in descending order to avoid index shifting
        const sorted = [...new Set(indices)].sort((a, b) => b - a);
        for (const idx of sorted) {
            if (idx >= 0 && idx < this.equations.length) {
                this.equations.splice(idx, 1);
            }
        }
    },

    getAllVariables() {
        const variables = new Set();
        for (const eq of this.equations) {
            Object.keys(eq.terms).forEach(v => variables.add(v));
        }
        return Array.from(variables).sort();
    },

    buildAugmentedMatrix() {
        const variables = this.getAllVariables();
        const matrix = [];

        for (const eq of this.equations) {
            const row = [];
            for (const v of variables) {
                row.push(eq.terms[v] || 0);
            }
            row.push(eq.constant);
            matrix.push(row);
        }

        return { matrix, variables };
    },

    gaussianElimination(matrix) {
        const m = matrix.length;
        const n = matrix[0].length;
        const augMatrix = matrix.map(row => [...row]);

        let pivotRow = 0;

        for (let col = 0; col < n - 1 && pivotRow < m; col++) {
            // Find pivot
            let maxRow = pivotRow;
            for (let row = pivotRow + 1; row < m; row++) {
                if (Math.abs(augMatrix[row][col]) > Math.abs(augMatrix[maxRow][col])) {
                    maxRow = row;
                }
            }

            if (Math.abs(augMatrix[maxRow][col]) < 1e-10) {
                continue; // Skip column if all zeros
            }

            // Swap rows
            [augMatrix[pivotRow], augMatrix[maxRow]] = [augMatrix[maxRow], augMatrix[pivotRow]];

            // Make pivot 1
            const pivot = augMatrix[pivotRow][col];
            for (let j = 0; j < n; j++) {
                augMatrix[pivotRow][j] /= pivot;
            }

            // Eliminate column
            for (let row = 0; row < m; row++) {
                if (row !== pivotRow) {
                    const factor = augMatrix[row][col];
                    for (let j = 0; j < n; j++) {
                        augMatrix[row][j] -= factor * augMatrix[pivotRow][j];
                    }
                }
            }

            pivotRow++;
        }

        return augMatrix;
    },

    solve() {
        if (this.equations.length === 0) {
            return { status: 'error', message: 'No equations to solve' };
        }

        const { matrix, variables } = this.buildAugmentedMatrix();
        const rref = this.gaussianElimination(matrix);
        
        const numVars = variables.length;
        const numEqs = rref.length;

        // Check for inconsistency
        for (let i = 0; i < numEqs; i++) {
            const allZero = rref[i].slice(0, numVars).every(x => Math.abs(x) < 1e-10);
            if (allZero && Math.abs(rref[i][numVars]) > 1e-10) {
                return { status: 'unsolvable', message: 'System has no solution (inconsistent)' };
            }
        }

        // Extract solution
        const solution = {};
        const free = [];
        const relations = [];

        for (let i = 0; i < numVars; i++) {
            const v = variables[i];
            let found = false;

            for (let row = 0; row < numEqs; row++) {
                if (Math.abs(rref[row][i]) > 1e-10) {
                    // Check if this is a leading 1
                    let isLeading = true;
                    for (let j = 0; j < i; j++) {
                        if (Math.abs(rref[row][j]) > 1e-10) {
                            isLeading = false;
                            break;
                        }
                    }

                    if (isLeading && Math.abs(rref[row][i] - 1) < 1e-10) {
                        // This variable has a leading 1
                        let value = rref[row][numVars];
                        const deps = [];

                        for (let j = i + 1; j < numVars; j++) {
                            if (Math.abs(rref[row][j]) > 1e-10) {
                                deps.push({ var: variables[j], coef: -rref[row][j] });
                            }
                        }

                        if (deps.length === 0) {
                            // Unique solution for this variable
                            solution[v] = this.applyConstraints(value);
                        } else {
                            // Parametric solution
                            const terms = deps.map(d => {
                                const coef = d.coef === 1 ? '' : (d.coef === -1 ? '-' : d.coef.toFixed(4));
                                return `${d.coef >= 0 ? '+' : ''}${coef}${d.var}`;
                            }).join('');
                            relations.push(`${v} = ${value.toFixed(4)}${terms}`);
                        }
                        found = true;
                        break;
                    }
                }
            }

            if (!found) {
                free.push(v);
            }
        }

        if (free.length > 0) {
            return {
                status: 'parametric',
                solution,
                relations,
                freeVariables: free,
                message: 'System has infinitely many solutions (underdetermined)'
            };
        }

        if (Object.keys(solution).length === variables.length) {
            return { status: 'unique', solution };
        }

        return { status: 'partial', solution, relations };
    },

    applyConstraints(value) {
        let result = value;

        if (this.flags.integer) {
            result = Math.round(result);
        }

        if (this.flags.positive && result < 0) {
            return null; // No valid solution under constraints
        }

        return result;
    },

    setFlag(flag, value) {
        if (flag in this.flags) {
            this.flags[flag] = value;
            return true;
        }
        return false;
    }
};
