// Implementation of the Beta probability density function
// Without the use of factorials to overcome the problem of getting inf values
// No external libraries required
// For naive implementation, see beta_naive.js
// Roy Hung 2019
// https://github.com/royhzq/betajs

export default function betaPDF(x: number, a: number, b: number) {
    // Beta probability density function impementation
    // using logarithms, no factorials involved.
    // Overcomes the problem with large integers
    return Math.exp(lnBetaPDF(x, a, b))
}

function lnBetaPDF(x: number, a: number, b: number) {
        // Log of the Beta Probability Density Function
    return ((a-1)*Math.log(x) + (b-1)*Math.log(1-x)) - lnBetaFunc(a,b)
}

function lnBetaFunc( a: number, b: number ) {
	// Log Beta Function
	// ln(Beta(x,y))
    let foo = 0.0;

    for( let i=0; i<a-2; i++ ) {
        foo += Math.log(a-1-i);
    }
    for( let i=0; i<b-2; i++ ) {
        foo += Math.log(b-1-i);
    }
    for( let i=0; i<a+b-2; i++ ) {
        foo -= Math.log(a+b-1-i);
    }
    return foo
}

/*function betaFunc( x: number, y: number ) {
    // Beta Function
    // Beta(x,y) = e^(ln(Beta(x,y))
    return Math.exp(lnBetaFunc(x,y));
}
*/