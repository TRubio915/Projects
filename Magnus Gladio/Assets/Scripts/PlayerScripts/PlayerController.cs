using System.Collections;
using System.Collections.Generic;
using System;
using UnityEngine;
using UnityEngine.SceneManagement;

public class PlayerController : MonoBehaviour
{
    private Rigidbody2D rb;

    private Animator anim;
    private MeleeAttackController m_AttackController;
    private RangedAttackWithTeleportation r_AttackController;
    private Material playerMat;

    [SerializeField]
    private GameObject cameraObject;
    private GameManager manager;


    private bool
      isFacingRight = true,
      isWalking,
      isRunning,
      isGrounded,
      isCloseToWall,
      isTouchedEnemy;
    public bool
        knockback;

    [Header("Movement Variables")]
    public bool canMove, canFlip;
    public bool isDead;
    public float
        movementSpeed;
    [SerializeField]
    private float
        maxSpeedWalking,
        maxSpeedRunning,
        accelerationFactor;

    private float upButtonDownTimer = 0;

    private float
        movementInputDirection, camPeekInputDirection;

    private int
        facingDirection;

    public Vector2 previousVelocity;

    [Header("Ledge Check")]
    public Transform ledgeCheck;

    private bool isTouchingLedge, canClimbLedge = false,
        ledgeDetected;

    
    [SerializeField]
    private float ledgeClimbxOffset1,
        ledgeClimbyOffset1,
        ledgeClimbxOffset2,
        ledgeClimbyOffset2;
    private Vector2 ledgePos1, ledgePos2, ledgePositionBot;

    [Header("Slippery Movement")]

    public LayerMask whatIsSlippery;
    public float slidingVelocityFactor;
    public float currentMaxSlidingSpeed;
    public float maxSlidingSpeed;
    private bool isSliding = false;

    [Space(10)]
    [Header("Smoke Particle")]
    [SerializeField]
    private float smokeSpawnTimeSet;
    private float smokeSpawnTime;
    [SerializeField]
    private Transform smokeSpawnPoint;
    [SerializeField]
    private GameObject smokeParticlePrefab;

    [Space(10)]
    [Header("Ground Check")]
    [SerializeField]
    private float
        groundCheckRadius;
    [SerializeField]
    private Transform
       groundCheck;
    [SerializeField]
    private LayerMask whatIsGround;


    [Space(10)]
    [Header("Wall Check")]
    [SerializeField]
    private float wallCheckRadius;
    [SerializeField]
    private Transform wallCheck;
    [SerializeField]
    private LayerMask whatIsWall;


    [Space(10)]
    [Header("KnockBack")]
    [SerializeField]
    private float knockBackDuration;
    [SerializeField]
    private Vector2 knockBackSpeed;
    private float
        knockBackStartTime;

    [Space(10)]
    [Header("Dash")]
    public bool dashUnlocked;
    public bool dashUpgraded;
    private bool
        isDashing = false,
        isUpgradedDash = false,
        canDash = true;

    [SerializeField]
    public float
        dashTime,
        dashSpeed,
        distanceBetweenImages,
        dashCoolDown,
        upgradedDashStaminaCost;
    private float
        dashTimeLeft,
        lastImageXPos,
        lastDash = -100;

    [Space(10)]
    [Header("WallSlide")]
    public float wallSlidingSpeed;
    public float baseWallSlidingSpeed;
    private bool isTouchingWall;

    [SerializeField]
    private float
        wallCheckDistance;

    public bool wallThrowUnlocked, isWallSliding;

    [Space(10)]
    [Header("Spell")]
    public GameObject spellPrefab;

    public bool spellAbilityUnlocked;
    private bool isUsingSpell = false;
    private float spellCastCD;
    public float spellCastRate;
    public float staminaCostForSpell;


    [Space(10)]
    [Header("Heal")]
    public float healRate;
    private bool isHealing = false;
    private float healCD;
    public float staminaCostForHeal;
    public float amountToHeal;
    private float buttonDownTimer;
    public float timeTakenToHeal;


    [Space(10)]
    [Header("Save/Load")]
    public LayerMask whatIsSavePoint;
    public LayerMask whatIsTeleporter;

    // Start is called before the first frame update
    void Start()
    {
        rb = GetComponent<Rigidbody2D>();
        anim = GetComponent<Animator>();
        m_AttackController = GetComponent<MeleeAttackController>();
        r_AttackController = GetComponent<RangedAttackWithTeleportation>();
        if (GameObject.FindGameObjectWithTag("Game Manager"))
            manager = GameObject.FindGameObjectWithTag("Game Manager").GetComponent<GameManager>();
        playerMat = GetComponent<SpriteRenderer>().material;
        facingDirection = 1;

        canMove = true;
        canFlip = true;
        isWalking = true;
        isRunning = false;
        isGrounded = true;
        isCloseToWall = false;
        baseWallSlidingSpeed = wallSlidingSpeed;

        currentMaxSlidingSpeed = maxSlidingSpeed;
    }

    // Update is called once per frame
    void Update()
    {

        if (!GameSettings.isGamePaused && !isDead && manager.isGameStart)
        {
            CheckInput();
            CheckMovementDirection();
            CheckSurroundings();
            CheckKnockBack();
            CheckAnimation();
            checkIfWallSliding();
            if (dashUnlocked)
                CheckDash();
            CheckLedgeClimb();
            CheckIfCanUseSpell();
            CheckIfCanUseHeal();
        }

    }
    private void FixedUpdate()
    {
        if (!isDead && manager.isGameStart)
        {
            if (canMove && !knockback)
            {
                ApplyMovement();
            }


            if (movementInputDirection == 0 && isGrounded && r_AttackController.isTeleporting != true && !m_AttackController.isUpAttacking)
            {
                upButtonDownTimer += Time.deltaTime;
                
            }
            else
            {
                upButtonDownTimer = 0;
                cameraObject.GetComponent<CameraController>().camPeek = false;
            }

            if(upButtonDownTimer >= 1)
            {
                cameraObject.GetComponent<CameraController>().camPeek = true;
                cameraObject.GetComponent<CameraController>().CamPeek(camPeekInputDirection);
            }
        }

    }


    private void CheckInput()
    {
        //Get Horizontal axis for movement Direction
        movementInputDirection = Input.GetAxisRaw("Horizontal");
        camPeekInputDirection = Input.GetAxis("Vertical");
        
        if (isWallSliding)
        {
            if (isFacingRight)
            {
                movementInputDirection = Mathf.Clamp(movementInputDirection, -1, 0);
            }
            else
            {
                movementInputDirection = Mathf.Clamp(movementInputDirection, 0, 1);
            }
        }
        else
        {
            movementInputDirection = Mathf.Clamp(movementInputDirection, -1, 1);
        }

        //Presses Shift Key to run
        if (Input.GetButtonDown("WalkRun") && isGrounded)
        {
            isWalking = false;
            isRunning = true;
        }
        //Release key to walk
        else if (Input.GetButtonUp("WalkRun") && isGrounded)
        {
            isWalking = true;
            isRunning = false;
        }

        if (Input.GetButtonDown("Dash") && canDash && dashUnlocked && !r_AttackController.isBeganThrowing)
        {
            if (Time.time >= lastDash + dashCoolDown)
            {
                if (dashUpgraded && PlayerStats.PlayerStamina > upgradedDashStaminaCost)
                {
                    AttemptToDash();
                }
                else if (!dashUpgraded)
                {
                    AttemptToDash();
                }
            }
        }



        if (Input.GetKeyDown(KeyCode.Alpha1))
        {            
            PlayerStats.PlayerMaxHealth = 1000;
            PlayerStats.PlayerHealth = 1000;
            PlayerStats.PlayerMaxStamina = 1000;
            PlayerStats.PlayerStamina = 1000;
        }
        if (Input.GetKeyDown(KeyCode.Alpha2))
        {
            manager.GoToNextArea();
        }
        if (Input.GetKeyDown(KeyCode.Alpha3))
        {
            manager.GoToNextRegion();
        }
        if (Input.GetKeyDown(KeyCode.Alpha4))
        {
            wallThrowUnlocked = true;
            dashUnlocked = true;
            dashUpgraded = true;
            spellAbilityUnlocked = true;            
        }

        if (Input.GetButtonDown("UseSpell") && spellAbilityUnlocked && !isUsingSpell && PlayerStats.PlayerStamina >= staminaCostForSpell)
        {
            UseSpell();
        }

        if (Input.GetButtonDown("Heal") && !isHealing && PlayerStats.PlayerStamina >= staminaCostForHeal && PlayerStats.PlayerHealth < PlayerStats.PlayerMaxHealth)
        {
            buttonDownTimer = Time.time;
            StartCoroutine(healEffect());

        }
        else if (Input.GetButton("Heal") && !isHealing && PlayerStats.PlayerStamina >= staminaCostForHeal && PlayerStats.PlayerHealth < PlayerStats.PlayerMaxHealth)
        {
            if (Time.time - buttonDownTimer > timeTakenToHeal)
            {
                buttonDownTimer = float.PositiveInfinity;
                playerMat.SetFloat("_Fade", 1);
                Heal();

            }
        }
        else
        {           
            buttonDownTimer = float.PositiveInfinity;
        }

        if(Input.GetButtonUp("Heal"))
        {
            StopAllCoroutines();
            playerMat.SetFloat("_Fade", 1);
        }

        if(Input.GetButtonDown("Interact") )
        {
            if (Physics2D.OverlapCircle(transform.position, 2, whatIsTeleporter))
            {
                manager.GoToNextArea();
            }



            if (Physics2D.OverlapCircle(transform.position, 2, whatIsSavePoint))
            {
                manager.OpenSaveMenu();
            }
        }


    }

    private void CheckMovementDirection()
    {

        if (isFacingRight && movementInputDirection < 0 && !isSliding)
        {
            Flip();
        }
        else if (!isFacingRight && movementInputDirection > 0 && !isSliding)
        {
            Flip();
        }

        if (Mathf.Abs(rb.velocity.x) >= 0.01f && isGrounded && !isWallSliding && !isSliding)
        {
            Sound[] sounds = FindObjectOfType<AudioManager>().sounds;
            Sound s = Array.Find(sounds, sound => sound.name == "Run");
            if (s.source.isPlaying == false)
                FindObjectOfType<AudioManager>().Play("Run");
        }

    }

    private void ApplyMovement()
    {

        if (!isSliding)
        {
            // If Idle
            isSliding = false;
            if (movementInputDirection == 0)
                movementSpeed = 0;
            else
            {
                if (Time.time >= smokeSpawnTime + smokeSpawnTimeSet && isGrounded && !isWallSliding)
                {
                    Instantiate(smokeParticlePrefab, smokeSpawnPoint.position, Quaternion.identity);
                    smokeSpawnTime = Time.time;
                }
                else
                {
                    smokeSpawnTime -= Time.fixedDeltaTime / Time.timeScale;
                }

            }
            //Increase movespeed based on the acceleration factor
            movementSpeed += accelerationFactor * 5 * Time.fixedUnscaledDeltaTime / Time.timeScale;

            if (isRunning)
            {
                //Clamp max possible speed to maxSpeedRunning while running
                movementSpeed = Mathf.Clamp(movementSpeed, 0, maxSpeedRunning / Time.timeScale);


            }
            else if (isWalking)
            {
                //Clamp max possible speed to maxSpeedWalking while walking
                movementSpeed = Mathf.Clamp(movementSpeed, 0, maxSpeedWalking / Time.timeScale);
            }

            //Change velocity of rb based on the input direction and the movementspeed
            rb.velocity = new Vector2(movementInputDirection * movementSpeed, rb.velocity.y);

            if (isWallSliding)
            {
                if (rb.velocity.y < -wallSlidingSpeed)
                {
                    rb.velocity = new Vector2(rb.velocity.x, -wallSlidingSpeed);
                }
            }


        }
        else
        {
            // Sliding

           
            float x = rb.velocity.x;
            x += Time.fixedUnscaledDeltaTime * 10 * facingDirection;
            rb.velocity = new Vector2( x , rb.velocity.y);
           
            Mathf.Clamp(rb.velocity.x, -maxSlidingSpeed, maxSlidingSpeed);
            
            Debug.Log(rb.velocity.x);
        }
    }

    private void CheckSurroundings()
    {
        isGrounded = Physics2D.OverlapCircle(groundCheck.position, groundCheckRadius, whatIsGround);


        if (Physics2D.OverlapCircle(groundCheck.position, groundCheckRadius, whatIsSlippery))
        {
            isGrounded = true;
            isSliding = true;
            rb.velocity = new Vector2(10 * facingDirection, rb.velocity.y);
            GetComponent<BoxCollider2D>().enabled = false;                      
        }
        else
        {
            isSliding = false;
            GetComponent<BoxCollider2D>().enabled = true;
        }

        if (!isGrounded)
        {
            isTouchingWall = Physics2D.Raycast(wallCheck.position, transform.right, wallCheckDistance, whatIsWall);
            GetComponent<BoxCollider2D>().enabled = true;
            isSliding = false;
        }
        else
            isTouchingWall = false;

        isCloseToWall = Physics2D.OverlapCircle(wallCheck.position, wallCheckRadius, whatIsWall);

        isTouchingLedge = Physics2D.Raycast(ledgeCheck.position, transform.right, wallCheckDistance, whatIsGround);

        if (isTouchingWall && !isTouchingLedge && !ledgeDetected && !isGrounded)
        {
            ledgeDetected = true;
            ledgePositionBot = wallCheck.position;
            canMove = false;
            canFlip = false;
            canDash = false;
        }
    }

    public void StopAllBasicMovement()
    {
        canMove = false;
        canFlip = false;

    }
    public void SetAllBasicMovement()
    {
        canMove = true;
        canFlip = true;

    }

    public void KnockBack(int direction)
    {
        knockback = true;
        canFlip = false;
        knockBackStartTime = Time.time;
        rb.velocity = new Vector2(knockBackSpeed.x * direction, knockBackSpeed.y);
    }

    private void CheckKnockBack()
    {
        if (Time.time >= knockBackStartTime + knockBackDuration && knockback)
        {
            //anim.SetBool("isHurt", false);
            knockback = false;
            canFlip = true;
            SetAllBasicMovement();
            rb.velocity = new Vector2(0.0f, rb.velocity.y);
        }
    }

    private void AttemptToDash()
    {
        isDashing = true;
        dashTimeLeft = dashTime;
        lastDash = Time.time;
        FindObjectOfType<AudioManager>().Play("Dash");
        PlayerAfterImagePool.Instance.GetFromPool();
        lastImageXPos = transform.position.x;

        if (dashUpgraded)
        {
            isUpgradedDash = true;
            manager.ReducePlayerStamina(upgradedDashStaminaCost);
        }
    }

    private void CheckDash()
    {
        if (isDashing)
        {
            if (dashTimeLeft > 0)
            {
                if (dashUpgraded)
                {
                    this.gameObject.layer = 14;                   
                }
                canMove = false;
                canFlip = false;
                //canjump = false;
                rb.velocity = new Vector2(facingDirection, 0) * dashSpeed * (1f / Time.timeScale);
                dashTimeLeft -= Time.deltaTime / Time.timeScale;

                if (Mathf.Abs(transform.position.x - lastImageXPos) > distanceBetweenImages)
                {
                    PlayerAfterImagePool.Instance.GetFromPool();
                    lastImageXPos = transform.position.x;
                }
            }

            if (dashTimeLeft <= 0 || isTouchedEnemy || isCloseToWall)
            {
                if (dashUpgraded)
                {
                    isUpgradedDash = false;
                    this.gameObject.layer = 8;
                    //GetComponent<BoxCollider2D>().enabled = true;
                    //GetComponent<CircleCollider2D>().enabled = true;
                }
                canMove = true;
                canFlip = true;
                isDashing = false;
                rb.velocity = Vector2.zero;
            }
        }
    }

    public bool GetDashStatus()
    {
        return isDashing;
    }
    public bool getUpgradedDashStatus()
    {
        return isUpgradedDash;
    }
    
    private void checkIfWallSliding()
    {
        if (isTouchingWall)
        {
            isWallSliding = true;
            //isGrounded = true;
        }
        else
        {
            isWallSliding = false;
        }
    }


    private void CheckIfCanUseSpell()
    {
        if (isUsingSpell && spellCastCD > 0)
        {
            spellCastCD -= Time.deltaTime;
        }
        if (spellCastCD <= 0)
        {
            isUsingSpell = false;
        }
    }
    private void UseSpell()
    {
        isUsingSpell = true;
        spellCastCD = spellCastRate;

        //Anim
        manager.ReducePlayerStamina(staminaCostForSpell);
        Instantiate(spellPrefab, transform.position, Quaternion.identity);

        anim.SetTrigger("CastSpell");

    }


    private void Heal()
    {
        if (PlayerStats.PlayerHealth < PlayerStats.PlayerMaxHealth)
        {
            isHealing = true;

            healCD = healRate;

            manager.ReducePlayerStamina(staminaCostForHeal);

            manager.IncreasePlayerhealth(amountToHeal);

        }

    }



    private void CheckIfCanUseHeal()
    {
        
        if (isHealing && healCD > 0)
        {
            healCD -= Time.deltaTime;
            
        }

        if(healCD <= 0)
        {
            isHealing = false;
            
        }
    }

    IEnumerator healEffect()
    {
        float fade = 0.25f;
        float t = 0;
        while (fade < 1)
        {
            yield return null;
            t += Time.deltaTime / timeTakenToHeal;
            
            fade = Mathf.Lerp(0.25f, 1, t);

           // fade += Time.deltaTime/3;
            playerMat.SetFloat("_Fade", fade);
           
        }

    }

    private void CheckLedgeClimb()
    {
        if (ledgeDetected && !canClimbLedge)
        {
            canClimbLedge = true;

            if (isFacingRight)
            {
                ledgePos1 = new Vector2(Mathf.Floor(ledgePositionBot.x + wallCheckDistance) - ledgeClimbxOffset1, Mathf.Floor(ledgePositionBot.y) + ledgeClimbyOffset1);
                ledgePos2 = new Vector2(Mathf.Floor(ledgePositionBot.x + wallCheckDistance) + ledgeClimbxOffset2, Mathf.Floor(ledgePositionBot.y) + ledgeClimbyOffset2);
            }
            else
            {
                ledgePos1 = new Vector2(Mathf.Ceil(ledgePositionBot.x - wallCheckDistance - 0.2f) + ledgeClimbxOffset1, Mathf.Floor(ledgePositionBot.y) + ledgeClimbyOffset1);
                ledgePos2 = new Vector2(Mathf.Ceil(ledgePositionBot.x - wallCheckDistance) - ledgeClimbxOffset2, Mathf.Floor(ledgePositionBot.y) + ledgeClimbyOffset2);
            }

            anim.SetBool("canClimbLedge", canClimbLedge);
            //FindObjectOfType<AudioManager>().Play("LedgeClimb");
            StartCoroutine("BeginLedgeClimb");
            StartCoroutine(FinishLedgeClimb());


        }

        if (canClimbLedge)
        {
            transform.position = ledgePos1;
            StopAllBasicMovement();           
        }
    }

    IEnumerator BeginLedgeClimb()
    {
        yield return new WaitForSeconds(0.5f);
        float t = 0;

        while( t < 0.5)
        {
            yield return null;
            transform.position = Vector2.Lerp(transform.position, new Vector2(transform.position.x, transform.position.y + 2), 0.5f);
            t += Time.deltaTime;
        }
       
    }
    IEnumerator FinishLedgeClimb()
    {
        yield return new WaitForSeconds(1f);
        canClimbLedge = false;
        //spriteObject.transform.localPosition = Vector3.zero;
        transform.position = ledgePos2;
        SetAllBasicMovement();
        ledgeDetected = false;
        anim.SetBool("canClimbLedge", canClimbLedge);
        
    }


    private void CheckAnimation()
    {
        anim.SetFloat("moveSpeed", rb.velocity.magnitude);
        anim.SetBool("isWallSliding", isWallSliding);
        anim.SetBool("isGrounded", isGrounded);
        anim.SetFloat("yVel", rb.velocity.y * -1);
        anim.SetBool("isSliding", isSliding);
    }

    public void Flip()
    {
        if (canFlip && !knockback)
        {

            isFacingRight = !isFacingRight;
            transform.Rotate(0, 180, 0f, 0);
            facingDirection *= -1;
        }

    }
    public int GetMovementDirectionOfPlayer()
    {
        return facingDirection;
    }
    public bool IsPlayerGrounded()
    {
        return isGrounded;
    }

    public bool isPlayerCloseToWall()
    {
        return isCloseToWall;
    }
    public bool isPlayerRunning()
    {
        return isRunning;
    }

    private void Damage(float[] attackDetails)
    {
        m_AttackController.Damage(attackDetails);
    }

    public void PlayerDied()
    {
        isDead = true;        
        anim.SetBool("isDead", true);
    }

    private void LoadGameOverScreen()
    {
        Destroy(manager);
        SceneManager.LoadScene(0);
    }

    private void OnDrawGizmos()
    {
        if (groundCheck != null)
        {
            Gizmos.DrawWireSphere(groundCheck.position, groundCheckRadius);
            Gizmos.color = Color.red;
            //Gizmos.DrawWireSphere(wallCheck.position, wallCheckRadius);
            Gizmos.DrawLine(wallCheck.position, new Vector2(wallCheck.position.x + wallCheckDistance, wallCheck.position.y));

            Gizmos.DrawLine(ledgeCheck.position, new Vector2(ledgeCheck.position.x + wallCheckDistance, ledgeCheck.position.y));
        }
    }


    private void OnCollisionEnter2D(Collision2D collision)
    {
        if (collision.gameObject.CompareTag("RegionCameraSwitcher"))
        {
            r_AttackController.resetThrowingVars();
            ScenePartLoader loader = GameObject.Find(manager.GetSceneName()).GetComponent<ScenePartLoader>();
        }
        if (collision.gameObject.CompareTag("Enemy"))
        {
            isTouchedEnemy = true;
        }
    }
    private void OnCollisionExit2D(Collision2D collision)
    {
        if (collision.gameObject.CompareTag("Enemy"))
        {
            isTouchedEnemy = false;
        }
    }

    private void OnTriggerEnter2D(Collider2D collision)
    {
        if (collision.CompareTag("SavePoint"))
        {
            //manager.SaveGame(GameSettings.saveSlot);
        }
    }
    public void speedUpgrade()
    {
        maxSpeedRunning += 2;
    }

}
