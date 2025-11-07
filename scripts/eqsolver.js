// Equation Solver Module
const EqSolver = {
    equations: [],
    equationStrings: [], // Store original equation strings
    isLinear: true,
    flags: {
        complex: false,
        integer: false,
        positive: false
    },

    reset() {
        this.equations = [];
        this.equationStrings = [];
        this.isLinear = true;
        this.flags = {
            complex: false,
            integer: false,
            positive: false,
            pi: false,
            fraction: false
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
        // Store original equation string
        this.equationStrings.push(eqStr);
        
        // Check if equation contains non-linear elements
        const nonlinearPatterns = [
            /\b(sin|cos|tan|exp|ln|log|sqrt|abs)\s*\(/i,  // Trig/math functions
            /\^/,                                           // Power operator
            /\*\s*[a-z]/i,                                 // Multiplication with variable
            /[a-z]\s*\*/i,                                 // Variable multiplication
            /[a-z]\s*[a-z]/i                              // Adjacent variables (like xy)
        ];
        
        const isNonlinear = nonlinearPatterns.some(pattern => pattern.test(eqStr));
        
        if (isNonlinear) {
            // Non-linear equation
            this.isLinear = false;
            this.equations.push({ type: 'nonlinear', expression: eqStr });
        } else {
            // Try to parse as linear equation
            try {
                const parsed = this.parseEquation(eqStr);
                this.equations.push(parsed);
            } catch (e) {
                // Failed to parse as linear, treat as non-linear
                this.isLinear = false;
                this.equations.push({ type: 'nonlinear', expression: eqStr });
            }
        }
        
        return this.equations.length - 1;
    },

    showEquations() {
        return this.equationStrings.map((eq, idx) => `[${idx}] ${eq}`);
    },

    deleteEquations(indices) {
        // Sort indices in descending order to avoid index shifting
        const sorted = [...new Set(indices)].sort((a, b) => b - a);
        for (const idx of sorted) {
            if (idx >= 0 && idx < this.equations.length) {
                this.equations.splice(idx, 1);
                this.equationStrings.splice(idx, 1);
            }
        }
        
        // Re-check if still linear
        this.isLinear = this.equations.every(eq => eq.type !== 'nonlinear');
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

        // Check if all equations are linear
        if (this.isLinear && this.equations.every(eq => eq.type !== 'nonlinear')) {
            return this.solveLinear();
        } else {
            return this.solveNonlinear();
        }
    },

    solveLinear() {
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
                            solution[v] = this.formatOutput(value);
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

    formatOutput(value) {
        if (value === null) return null;

        // Apply constraints first
        value = this.applyConstraints(value);
        if (value === null) return null;

        // Format as PI if flag is set
        if (this.flags.pi) {
            const piRatio = value / Math.PI;
            
            // Check if it's close to a simple fraction of PI
            const tolerance = 1e-6;
            const denominators = [1, 2, 3, 4, 6, 8, 12];
            
            for (const denom of denominators) {
                const numerator = Math.round(piRatio * denom);
                if (Math.abs(piRatio - numerator / denom) < tolerance) {
                    if (numerator === 0) return '0';
                    if (numerator === denom) return 'π';
                    if (numerator === -denom) return '-π';
                    if (denom === 1) return `${numerator}π`;
                    
                    const sign = numerator < 0 ? '-' : '';
                    const absNum = Math.abs(numerator);
                    if (absNum === 1) return `${sign}π/${denom}`;
                    return `${sign}${absNum}π/${denom}`;
                }
            }
            
            // Not a simple fraction of PI, show decimal with π notation
            return `${piRatio.toFixed(6)}π`;
        }

        // Format as fraction if flag is set
        if (this.flags.fraction) {
            return this.toFraction(value);
        }

        return value;
    },

    toFraction(decimal, maxDenominator = 10000) {
        if (Math.abs(decimal - Math.round(decimal)) < 1e-9) {
            return Math.round(decimal).toString();
        }

        const sign = decimal < 0 ? '-' : '';
        decimal = Math.abs(decimal);

        let bestNumerator = 1;
        let bestDenominator = 1;
        let minError = Math.abs(decimal - 1);

        for (let denominator = 1; denominator <= maxDenominator; denominator++) {
            const numerator = Math.round(decimal * denominator);
            const error = Math.abs(decimal - numerator / denominator);
            
            if (error < minError) {
                minError = error;
                bestNumerator = numerator;
                bestDenominator = denominator;
            }
            
            if (error < 1e-9) break;
        }

        if (bestDenominator === 1) {
            return `${sign}${bestNumerator}`;
        }

        return `${sign}${bestNumerator}/${bestDenominator}`;
    },

    // Nonlinear equation solver using Newton-Raphson method
    solveNonlinear() {
        try {
            // Extract all variables from equations
            const variables = this.extractVariables();
            
            if (variables.length === 0) {
                return { status: 'error', message: 'No variables found in equations' };
            }

            // For single equation, use Newton-Raphson
            if (this.equationStrings.length === 1 && variables.length === 1) {
                return this.solveNewtonRaphson(this.equationStrings[0], variables[0]);
            }

            // For multiple equations/variables, use multidimensional Newton
            return this.solveNewtonMulti(variables);
        } catch (error) {
            return { status: 'error', message: error.message };
        }
    },

    extractVariables() {
        const varSet = new Set();
        const varRegex = /\b([a-z])\b/gi;
        
        for (const eqStr of this.equationStrings) {
            const matches = eqStr.matchAll(varRegex);
            for (const match of matches) {
                const v = match[1].toLowerCase();
                // Exclude function names
                if (!['sin', 'cos', 'tan', 'log', 'ln', 'exp', 'abs', 'sqrt', 'e', 'pi'].includes(v)) {
                    varSet.add(v);
                }
            }
        }
        
        return Array.from(varSet).sort();
    },

    // Evaluate expression using Function constructor (safe for math expressions)
    evaluateExpression(expr, varValues) {
        try {
            // Replace variables with values
            let evalExpr = expr;
            for (const [varName, value] of Object.entries(varValues)) {
                // Use word boundaries to avoid replacing parts of function names
                const regex = new RegExp(`\\b${varName}\\b`, 'g');
                evalExpr = evalExpr.replace(regex, `(${value})`);
            }

            // Replace math functions with Math. equivalents
            evalExpr = evalExpr
                .replace(/\bsin\(/g, 'Math.sin(')
                .replace(/\bcos\(/g, 'Math.cos(')
                .replace(/\btan\(/g, 'Math.tan(')
                .replace(/\blog\(/g, 'Math.log10(')
                .replace(/\bln\(/g, 'Math.log(')
                .replace(/\bexp\(/g, 'Math.exp(')
                .replace(/\babs\(/g, 'Math.abs(')
                .replace(/\bsqrt\(/g, 'Math.sqrt(')
                .replace(/\bpi\b/g, 'Math.PI')
                .replace(/\be\b/g, 'Math.E')
                .replace(/\^/g, '**'); // Power operator

            // Evaluate
            return Function(`"use strict"; return (${evalExpr})`)();
        } catch (e) {
            throw new Error(`Failed to evaluate expression: ${e.message}`);
        }
    },

    // Numerical derivative using central difference
    derivative(expr, varName, varValues, h = 1e-6) {
        const values1 = { ...varValues, [varName]: varValues[varName] + h };
        const values2 = { ...varValues, [varName]: varValues[varName] - h };
        
        const f1 = this.evaluateExpression(expr, values1);
        const f2 = this.evaluateExpression(expr, values2);
        
        return (f1 - f2) / (2 * h);
    },

    // Single variable Newton-Raphson
    solveNewtonRaphson(equation, varName, maxIterations = 100, tolerance = 1e-6) {
        // Split equation into LHS and RHS
        const sides = equation.split('=');
        if (sides.length !== 2) {
            return { status: 'error', message: 'Invalid equation format' };
        }

        // Create f(x) = LHS - RHS
        const f = `(${sides[0].trim()}) - (${sides[1].trim()})`;

        // Initial guess
        let x = 0;
        let iterations = 0;

        while (iterations < maxIterations) {
            const varValues = { [varName]: x };
            
            // Evaluate f(x)
            const fx = this.evaluateExpression(f, varValues);
            
            if (Math.abs(fx) < tolerance) {
                // Found solution
                const result = this.formatOutput(x);
                if (result === null) {
                    return { status: 'error', message: 'Solution violates constraints' };
                }
                return {
                    status: 'unique',
                    solution: { [varName]: result },
                    iterations,
                    method: 'Newton-Raphson'
                };
            }

            // Evaluate f'(x)
            const fpx = this.derivative(f, varName, varValues);
            
            if (Math.abs(fpx) < 1e-12) {
                // Derivative too small, try different starting point
                x = Math.random() * 10 - 5;
                iterations++;
                continue;
            }

            // Newton update: x = x - f(x)/f'(x)
            x = x - fx / fpx;
            iterations++;
        }

        return {
            status: 'error',
            message: `No convergence after ${maxIterations} iterations. Try different initial conditions.`
        };
    },

    // Multi-variable Newton method
    solveNewtonMulti(variables, maxIterations = 100, tolerance = 1e-6) {
        const n = variables.length;
        const m = this.equationStrings.length;

        if (n !== m) {
            return {
                status: 'error',
                message: `System is ${n > m ? 'underdetermined' : 'overdetermined'} (${n} variables, ${m} equations)`
            };
        }

        // Initial guess (all zeros)
        const x = variables.reduce((obj, v) => ({ ...obj, [v]: 0 }), {});
        
        for (let iter = 0; iter < maxIterations; iter++) {
            // Build Jacobian matrix and function vector
            const J = [];
            const F = [];

            for (let i = 0; i < m; i++) {
                const eq = this.equationStrings[i];
                const sides = eq.split('=');
                const f = `(${sides[0].trim()}) - (${sides[1].trim()})`;

                // Evaluate f_i(x)
                const fi = this.evaluateExpression(f, x);
                F.push(fi);

                // Compute Jacobian row (derivatives with respect to each variable)
                const jRow = [];
                for (const v of variables) {
                    const deriv = this.derivative(f, v, x);
                    jRow.push(deriv);
                }
                J.push(jRow);
            }

            // Check convergence
            const norm = Math.sqrt(F.reduce((sum, fi) => sum + fi * fi, 0));
            if (norm < tolerance) {
                // Apply constraints and format
                const solution = {};
                for (const v of variables) {
                    const val = this.formatOutput(x[v]);
                    if (val === null) {
                        return { status: 'error', message: `Solution for ${v} violates constraints` };
                    }
                    solution[v] = val;
                }
                
                return {
                    status: 'unique',
                    solution,
                    iterations: iter,
                    method: 'Multidimensional Newton'
                };
            }

            // Solve J * delta = -F using Gaussian elimination
            const augmented = J.map((row, i) => [...row, -F[i]]);
            const delta = this.solveLinearSystem(augmented);

            if (delta === null) {
                return { status: 'error', message: 'Jacobian is singular, cannot continue' };
            }

            // Update x = x + delta
            for (let i = 0; i < n; i++) {
                x[variables[i]] += delta[i];
            }
        }

        return {
            status: 'error',
            message: `No convergence after ${maxIterations} iterations`
        };
    },

    // Solve linear system Ax = b (A is augmented matrix)
    solveLinearSystem(augmented) {
        const n = augmented.length;
        const rref = this.gaussianElimination(augmented);
        
        // Back substitution
        const x = new Array(n);
        for (let i = n - 1; i >= 0; i--) {
            if (Math.abs(rref[i][i]) < 1e-10) {
                return null; // Singular matrix
            }
            x[i] = rref[i][n];
        }
        
        return x;
    }
};
