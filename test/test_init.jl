function test_init()

    @test CTDirectShooting._init_vector(1, 0.1, 3) == [0.0, 1.0, 0.1]
    
end