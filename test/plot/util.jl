# Tests only date_to_float currently

dates = [quartertodate("2019-Q1"),quartertodate("2019-Q2"),quartertodate("2019-Q3"),quartertodate("2019-Q4")]

output = date_to_float(dates, 2)
answer = [2019.25 2019.5 2019.75 2020.;
          2019.25 2019.5 2019.75 2020.]
@testset "Check date_to_float converts dates correctly" begin
    @test output == answer
end
