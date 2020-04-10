$(document).on('shiny:busy', function() {
  var $inputs = $('button,input');
console.log($inputs);
$inputs.prop('disabled', true);
});

$(document).on('shiny:idle', function() {
var $inputs = $('button,input');
console.log($inputs);
$inputs.prop('disabled', false);
})
