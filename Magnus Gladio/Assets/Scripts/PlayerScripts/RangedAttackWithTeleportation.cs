using System.Collections;
using System.Collections.Generic;
using System;
using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.Rendering.Universal;
using UnityEngine.UI;

public class RangedAttackWithTeleportation : MonoBehaviour
{

    private PlayerController controller;
    private MeleeAttackController m_AttackController;
    private Rigidbody2D rb;
    private Animator anim;
    private GameManager manager;

    [SerializeField]
    private GameObject crossHair;

    private bool isWeaponInHand = true;
    private bool outOfRange = false;
    private bool enemyStun = false;
    public bool isBeganThrowing = false;
    private bool isReleasing = false;
    
    public bool canThrow, isTeleporting, retrieveWeapon = false;
    int cameraSwitch = -1, regionSwitch = -1;

    Camera cam;

    public float jumpForce;

    [Space(10)]
    [Header("Camera Zoom variables")]

    [SerializeField]
    public float minCamZoomSize;
    [SerializeField]
    public float maxCamZoomSize;

    [Space(10)]
    [Header("Projectile Launch Variables")]
    
    [SerializeField]
    private GameObject projectile;
    
    [SerializeField] private float launchRate;

    private float launchCoolDown;

    [SerializeField]
    private Transform projectileSpawnPoint;

    private Vector3 aimPos;
    
    public float firingHeight = 10;
    public float maxFiringHeight;
    public float minFiringHeight;
    public float arrowForce;
    public LayerMask whatIsEnd, whatIsEnemy;
    private GameObject enemyTargeted;
    public GameObject targetIndicator;       
    private Vector2 dir;
    private Vector3 targetPos;
    private List<Vector3> pointList;
    [SerializeField]
    private LineRenderer projectilePath;
    
    [Space(10)]
    [SerializeField]
    private GameObject smokePrefab;

    public Volume globalVolume;
    DepthOfField dofComp;
    public CameraShaker cameraShaker;


    [Space(10)]
    [SerializeField]
    private Image arcMode, straightMode;
    public int currentMode;

    private Coroutine zoomout;


    // Start is called before the first frame update
    void Start()
    {
        Cursor.visible = false;
        pointList = new List<Vector3>();
        controller = GetComponent<PlayerController>();
        m_AttackController = GetComponent<MeleeAttackController>();
        if(GameObject.FindGameObjectWithTag("Game Manager"))
            manager = GameObject.FindGameObjectWithTag("Game Manager").GetComponent<GameManager>();
        rb = GetComponent<Rigidbody2D>();
        anim = GetComponent<Animator>();
        cam = Camera.main;
        cam.orthographicSize = minCamZoomSize;
        Volume volume = globalVolume.GetComponent<Volume>();
        
        DepthOfField tmpdof;
        
        if (volume.profile.TryGet<DepthOfField>(out tmpdof))
        {
            dofComp = tmpdof;
        }
        currentMode = 0;
        arcMode.color = new Color(1, 1, 1, 1);
        straightMode.color = new Color(1, 1, 1, 35 / 255);

    }

    public void SetupRangeAttackController()
    {
        if (SaveSystem.getFileNames().Count <= 0)
        {
            canThrow = false;
        }
        else
        {
            canThrow = true;
        }
    }

    // Update is called once per frame
    void Update()
    {
        if (!GameSettings.isGamePaused && !controller.isDead)
        {
            
            aimPos = Camera.main.ScreenToWorldPoint(Input.mousePosition);
            aimPos.z = 0;
            crossHair.transform.position = aimPos;

            CheckInput();
            
            CheckMode();

            if (retrieveWeapon)
            {                
                RetrieveWeapon();
            }
         
            if(!isWeaponInHand)
                projectilePath.positionCount = 0;
        }
    }
    private void CheckInput()
    {
        if (Input.GetButtonDown("Fire2"))
        {
            if (!outOfRange && isWeaponInHand && !isBeganThrowing)
            {
                BeginThrowingWeapon();
            }
            else if(!isWeaponInHand && !retrieveWeapon)
            {
                BeginTeleport();
            }
            
        }

        if (Input.GetButtonUp("Fire2") && isBeganThrowing)
        {
            if (isWeaponInHand && canThrow)
            {
                if(controller.IsPlayerGrounded() && !controller.isPlayerCloseToWall())
                {
                    
                        isReleasing = true;
                        StopCoroutine("BeginCookingProjectile");                    
                        
                        anim.SetTrigger("isThrowing");
                    
                   
                }
                else if(controller.isWallSliding && controller.wallThrowUnlocked)
                {
                    //isBeganThrowing = false;
                    
                        isReleasing = true;
                        StopCoroutine("BeginCookingProjectile");
                        FindObjectOfType<AudioManager>().Play("Throw");
                        ReleaseProjectile();
                    
                }
                else if(!controller.IsPlayerGrounded() && !controller.isPlayerCloseToWall())
                {
                    isReleasing = true;
                    StopCoroutine("BeginCookingProjectile");
                    FindObjectOfType<AudioManager>().Play("Throw");
                    ReleaseProjectile();
                }
                
            }            
        }

        if (Input.GetButtonDown("RetrieveWeapon"))
        {               
            if (!isWeaponInHand && retrieveWeapon == false)
            {
                retrieveWeapon = true;                
            }
        }

        if(Input.GetButtonDown("SwitchMode") && manager.isKunaiModeUnlocked)
        {
            if (currentMode == 0)
            {
                currentMode = 1;
                straightMode.color = new Color(1,1,1,1);
                arcMode.color = new Color(1,1,1,35 / 255);
            }
            else if (currentMode == 1)
            {
                currentMode = 0;
                arcMode.color = new Color(1, 1, 1, 1);
                straightMode.color = new Color(1, 1, 1, 35 / 255);
            }
            else
                currentMode = 0;
        }


        if(Input.GetButtonDown("Jump"))
        {
            if (!isWeaponInHand || !canThrow || !controller.IsPlayerGrounded()) return;

            GameObject projectileClone = Instantiate(projectile, projectileSpawnPoint.transform.position, Quaternion.identity);
            
            canThrow = false;
            isWeaponInHand = false;
            launchCoolDown = launchRate;
            if(controller.GetComponent<Rigidbody2D>().velocity == Vector2.zero)
                projectileClone.GetComponent<Rigidbody2D>().AddForce(Vector2.up * jumpForce);
            else
            {
                controller.previousVelocity = rb.velocity;
                if (controller.GetMovementDirectionOfPlayer() > 0)
                {
                    
                    Vector2 dir = new Vector2(1, 1);
                    Quaternion rotation = Quaternion.Slerp(projectileClone.transform.rotation, Quaternion.LookRotation(Vector3.forward, dir), Time.fixedDeltaTime * 100);
                    projectileClone.transform.rotation = rotation;
                    projectileClone.GetComponent<Rigidbody2D>().AddForce(dir * jumpForce);
                }
                else
                {
                    Vector2 dir = new Vector2(-1, 1);
                    Quaternion rotation = Quaternion.Slerp(projectileClone.transform.rotation, Quaternion.LookRotation(Vector3.forward, dir), Time.fixedDeltaTime * 100);
                    projectileClone.transform.rotation = rotation;
                    projectileClone.GetComponent<Rigidbody2D>().AddForce(dir * jumpForce);
                }

            }


            Invoke("Jump", 0.25f);
        }
    }

    private void Jump()
    {
        BeginTeleport();
        Vector2 target = Vector2.zero;

        GameObject projectile = GameObject.FindGameObjectWithTag("Projectile");
        if (projectile != null)
            target = GameObject.FindGameObjectWithTag("Projectile").transform.position;

        while (Vector2.Distance(transform.position, target) > 0.25f && projectile != null)
        {
            transform.position = Vector2.MoveTowards(transform.position, target, Time.unscaledDeltaTime);            
        }
        StartCoroutine(resetIndicatorIcon());
        isWeaponInHand = true;
        canThrow = true;
        isBeganThrowing = false;
        anim.SetBool("isTeleporting", false);
        if (projectile != null)
            Destroy(GameObject.FindGameObjectWithTag("Projectile").gameObject);
        m_AttackController.canAttack = true;
        isTeleporting = false;
        if (cameraSwitch != -1)
        {
            manager.switchCameraBoundsBetweenRegion(cameraSwitch);
            cameraSwitch = -1;
        }
        if (regionSwitch != -1)
        {
            manager.setCameraBounds(regionSwitch);
            regionSwitch = -1;
        }
    }

    private void CheckMode()
    {
        if (!isBeganThrowing)
        {
            if (currentMode == 0)
            {
                RaycastHit2D hitDown = Physics2D.Raycast(crossHair.transform.position, Vector2.down, 15, whatIsEnd);

                if (hitDown)
                {
                    //targetIndicator.GetComponent<SpriteRenderer>().color = new Color(1, 1, 0, 0.44f);
                    targetIndicator.transform.position = hitDown.point;
                    outOfRange = false;
                    //canThrow = true;
                }
                else
                {
                    //targetIndicator.GetComponent<SpriteRenderer>().color = new Color(1, 0, 0, 0.44f);
                    targetIndicator.transform.position = crossHair.transform.position;
                    if (isWeaponInHand && isBeganThrowing && controller.IsPlayerGrounded() && !controller.isPlayerCloseToWall())
                    {
                        isBeganThrowing = false;
                        StopCoroutine("BeginCookingProjectile");
                        anim.SetTrigger("isThrowing");
                    }
                    outOfRange = true;
                    //canThrow = false;
                }
            }
            else
            {
                targetIndicator.transform.position = crossHair.transform.position;                    
            }
            /*
            if (!enemyStun)
            {
                Collider2D hitEnemy = Physics2D.OverlapPoint(crossHair.transform.position, whatIsEnemy);

                if (hitEnemy != null && hitEnemy.GetComponent<Enemy>() != null)
                {
                    hitEnemy.GetComponent<Enemy>().Targeted();
                    enemyStun = true;
                    enemyTargeted = hitEnemy.gameObject;
                }
                else
                {
                    if (enemyTargeted != null)
                    {
                        enemyTargeted.GetComponent<Enemy>().UnTarget();
                    }
                }
            }
            else
            {
                if (enemyTargeted != null)
                {
                    Collider2D hitEnemy = Physics2D.OverlapPoint(crossHair.transform.position, whatIsEnemy);
                    if (hitEnemy == null)
                    {
                        enemyTargeted.GetComponent<Enemy>().UnTarget();
                        enemyTargeted = null;
                        enemyStun = false;
                    }
                }
            }*/
        }
    }

    private void BeginThrowingWeapon()
    {
        
        if (canThrow && launchCoolDown <= 0 &&( ((controller.IsPlayerGrounded() && !controller.isPlayerCloseToWall()) || (controller.isWallSliding && controller.wallThrowUnlocked))))
        {
            isBeganThrowing = true;
            if (!controller.isWallSliding)
            {
                if (aimPos.x > transform.position.x && controller.GetMovementDirectionOfPlayer() < 0)
                {
                    controller.Flip();
                }
                else if (aimPos.x < transform.position.x && controller.GetMovementDirectionOfPlayer() > 0)
                {
                    controller.Flip();
                }
            }
            dir = (aimPos - gameObject.transform.position).normalized;

            //rb.velocity = Vector2.zero;
            //controller.StopAllBasicMovement();
            targetPos = targetIndicator.transform.position;

            if (controller.isWallSliding && controller.wallThrowUnlocked)
            {
                Physics2D.gravity = new Vector2(Physics2D.gravity.x, 0);
                controller.wallSlidingSpeed = 0;
            }
            else
            {
                Physics2D.gravity = new Vector2(Physics2D.gravity.x, -98);
            }
            Time.timeScale = 0.1f;
            Time.fixedDeltaTime = 0.02F * Time.timeScale;
            rb.gravityScale = 2 / Time.timeScale;
            m_AttackController.isAttacking = false;
           
            if (isBeganThrowing && currentMode == 0)
                StartCoroutine(BeginCookingProjectile());

            /*if (enemyStun)
                ReleaseProjectile();*/
        }
    
    
    }
    private void BeginTeleport()
    {
        if (isTeleporting == false)
        {
            //Teleport to weapon                   
            Physics2D.gravity = new Vector2(Physics2D.gravity.x, -9.8f);            
            isTeleporting = true;
            Sound[] sounds = FindObjectOfType<AudioManager>().sounds;
            Sound s = Array.Find(sounds, sound => sound.name == "Teleport");
            if (s.source.isPlaying == false)
                FindObjectOfType<AudioManager>().Play("Teleport");

            anim.SetBool("isTeleporting", true);
            
           // StartCoroutine(cameraShaker.Shake(0.1f, 0.05f));
        }        
    }

    private IEnumerator BeginCookingProjectile()
    {   
        
        while (Input.GetButton("Fire2") && isBeganThrowing && canThrow && !GameSettings.isGamePaused)
        {
            m_AttackController.isAttacking = false;
            cam.orthographicSize = Mathf.MoveTowards(cam.orthographicSize, maxCamZoomSize, 5f * Time.deltaTime / Time.timeScale);
            launchCoolDown = launchRate;
            firingHeight = Mathf.MoveTowards(firingHeight, maxFiringHeight, 7.5f * Time.deltaTime / Time.timeScale);
            firingHeight = Mathf.Clamp(firingHeight, minFiringHeight, maxFiringHeight);

            
            dofComp.mode.value = DepthOfFieldMode.Bokeh;

            yield return null;

            if(currentMode == 0)
                simulateArc(targetPos, 0.65f);
            
        }
        canThrow = false;
    }
    private void ReleaseProjectile()
    {
        if (isReleasing == true)
        {
            Time.timeScale = 1;
            Time.fixedDeltaTime = 0.02F * Time.timeScale;
            Physics2D.gravity = new Vector2(Physics2D.gravity.x, -9.8f);
            if (rb != null)
                rb.gravityScale = 2;

            canThrow = false;
            FindObjectOfType<AudioManager>().Play("Throw");
            //StartCoroutine(cameraShaker.Shake(0.1f, 0.05f));
            launchCoolDown = launchRate;
            isWeaponInHand = false;

            m_AttackController.canAttack = false;
            manager.throwableIndicator.fillAmount = 0;
            controller.SetAllBasicMovement();
            int spawnDirection = controller.GetMovementDirectionOfPlayer() * -1;
            Vector3 pos = projectileSpawnPoint.position;
            if (controller.isWallSliding && controller.wallThrowUnlocked)
                pos = new Vector3(projectileSpawnPoint.position.x + spawnDirection, projectileSpawnPoint.position.y, projectileSpawnPoint.position.z);
            GameObject projectileClone = Instantiate(projectile, pos, Quaternion.identity);
            float rot = Mathf.Atan2(dir.y, dir.x) * Mathf.Rad2Deg;
            projectileClone.transform.rotation = Quaternion.Euler(0, 0, rot - 90);
            projectilePath.positionCount = 0;
            if(currentMode == 0)
                StartCoroutine(simulateProjectile(targetPos, 0.65f));
            else
            {
                controller.currentMaxSlidingSpeed = 5;
                GameObject projClone = GameObject.FindGameObjectWithTag("Projectile");
                float angle;
                angle = Mathf.Atan2(targetPos.y, targetPos.x) * Mathf.Rad2Deg - 90;
                projClone.GetComponent<Rigidbody2D>().AddForce((targetPos - projClone.transform.position).normalized * arrowForce, ForceMode2D.Impulse);
                controller.previousVelocity = projClone.GetComponent<Rigidbody2D>().velocity;
                zoomout = StartCoroutine(ZoomOutCamera());
            }
            if(currentMode == 0)
                StartCoroutine(resetCamera());
            isReleasing = false;
        }
    }

    IEnumerator simulateProjectile(Vector3 dest, float time)
    {
        GameObject projClone = GameObject.FindGameObjectWithTag("Projectile");
        var startPos = projClone.transform.position;
        var timer = 0.0f;
        dest = new Vector3(dest.x, dest.y, 0);
        Vector2 vel = Vector2.zero;
        //if (currentMode == 0)
        //{
            while (projClone != null && timer <= 1 && projClone.GetComponent<Arrow>().hasHit != true && retrieveWeapon != true && isTeleporting == false)
            {
                var height = Mathf.Sin(Mathf.PI * timer) * firingHeight;
                Transform target;
                Vector3 nextPos;
                Vector3 thisPos;
                float angle;
                nextPos = Vector3.Lerp(startPos, dest, timer) + Vector3.up * height;
                thisPos = projClone.transform.position;
                nextPos.x = nextPos.x - thisPos.x;
                nextPos.y = nextPos.y - thisPos.y;
                Quaternion rotation = Quaternion.Slerp(projClone.transform.rotation, Quaternion.LookRotation(Vector3.forward,nextPos), Time.fixedDeltaTime * 100);
                projClone.transform.rotation = rotation;
                //angle = Mathf.Atan2(nextPos.y, nextPos.x) * Mathf.Rad2Deg - 90;
                //projClone.transform.rotation = Quaternion.Euler(new Vector3(0, 0, angle));
                
                projClone.GetComponent<Rigidbody2D>().position = Vector3.Lerp(startPos, dest, timer) + Vector3.up * height;
                timer += Time.unscaledDeltaTime / time;
                controller.previousVelocity = projClone.GetComponent<Rigidbody2D>().velocity;
                Collider2D collider = Physics2D.OverlapPoint(projClone.GetComponent<Rigidbody2D>().position, whatIsEnd);

                if (collider != null)
                {
                    break;
                }
                yield return null;
            }
        
       /* else
        {
                yield return null;
                float angle;               
                angle = Mathf.Atan2(dest.y, dest.x) * Mathf.Rad2Deg - 90;
                projClone.transform.rotation = Quaternion.Euler(new Vector3(0, 0, -angle));
                projClone.GetComponent<Rigidbody2D>().velocity = (dest - projClone.transform.position).normalized * arrowForce * Time.unscaledDeltaTime;             

        }*/
        launchCoolDown = launchRate;
    }

    private void simulateArc(Vector3 dest, float time )
    {
        var startPos = projectileSpawnPoint.transform.position;
        var timer = 0.0f;
        dest = new Vector3(dest.x, dest.y, 0);
        pointList.Clear();
        projectilePath.transform.position = projectileSpawnPoint.position;

        if(controller.isWallSliding)
        {
            startPos = new Vector3(projectileSpawnPoint.position.x + (controller.GetMovementDirectionOfPlayer() * -1),
                                                    projectileSpawnPoint.position.y, projectileSpawnPoint.position.z);
            projectilePath.transform.position = new Vector3(projectileSpawnPoint.position.x + (controller.GetMovementDirectionOfPlayer() * -1), 
                                                    projectileSpawnPoint.position.y, projectileSpawnPoint.position.z);
        }
        Vector3 newPos = Vector3.zero ;
        while (timer <= 1.0f)
        {
            var height = Mathf.Sin(Mathf.PI * timer) * firingHeight;
            newPos = Vector3.Lerp(startPos, dest, timer) + Vector3.up * height; 
           
            timer += Time.unscaledDeltaTime / time;
            pointList.Add(newPos);
            Collider2D collider = Physics2D.OverlapPoint(newPos, whatIsEnd);

            if(collider != null)
            {
                break;
            }
        }
        projectilePath.positionCount = pointList.Count;
        projectilePath.SetPositions(pointList.ToArray());        
    }


    private IEnumerator ZoomOutCamera()
    {
        while (cam.orthographicSize < maxCamZoomSize)
        {
            yield return null;

            cam.orthographicSize = Mathf.MoveTowards(cam.orthographicSize, maxCamZoomSize, 10f * Time.unscaledDeltaTime);
        }
    }

    private void PlayerTeleported()
    {   

        Time.timeScale = 1;
        if (zoomout != null)
            StopCoroutine(zoomout);
        StartCoroutine(resetCamera());
        Time.fixedDeltaTime = 0.02F * Time.timeScale;
        Physics2D.gravity = new Vector2(Physics2D.gravity.x, -9.8f);
        controller.wallSlidingSpeed = controller.baseWallSlidingSpeed;
        rb.gravityScale = 2 ;
        firingHeight = minFiringHeight;        
        isBeganThrowing = false;
        dofComp.mode.value = DepthOfFieldMode.Off;
        StartCoroutine(moveTowardsProjectile());

    }

    private void RetrieveWeapon()
    {
        dofComp.mode.value = DepthOfFieldMode.Off;
        controller.wallSlidingSpeed = 1f;
        Sound[] sounds = FindObjectOfType<AudioManager>().sounds;
        Sound s = Array.Find(sounds, sound => sound.name == "Retrieve");
        if (s.source.isPlaying == false)
            FindObjectOfType<AudioManager>().Play("Retrieve");
        Time.timeScale = 1;
        Time.fixedDeltaTime = 0.02F * Time.timeScale;
        Physics2D.gravity = new Vector2(Physics2D.gravity.x, -9.8f);
        rb.gravityScale = 2;
        firingHeight = minFiringHeight;
        isBeganThrowing = false;
        GameObject projClone = GameObject.FindGameObjectWithTag("Projectile");
        projClone.GetComponent<Rigidbody2D>().velocity = Vector2.zero;
        projClone.GetComponent<BoxCollider2D>().enabled = false;
        projClone.GetComponent<Arrow>().isCalledBack = true;

        Vector3 targ = projClone.transform.position;
        targ.z = 0f;
        Vector3 objectPos = transform.position;
        targ.x = objectPos.x - targ.x;
        targ.y = objectPos.y - targ.y;
        if (zoomout != null)
            StopCoroutine(zoomout);
        StartCoroutine(resetCamera());
        float angle = Mathf.Atan2(targ.y, targ.x) * Mathf.Rad2Deg - 90;

        projClone.transform.rotation = Quaternion.Euler(new Vector3(0, 0, angle));

        if (Vector2.Distance(projClone.transform.position, transform.position) > 0.5f)
            projClone.transform.position = Vector3.MoveTowards(projClone.transform.position, transform.position, 2.0f);
        else
        {
            retrieveWeapon = false;
            StartCoroutine(resetIndicatorIcon());
            isWeaponInHand = true;
            canThrow = true;
            Destroy(projClone);
            m_AttackController.canAttack = true;
        }
    }

    private IEnumerator moveTowardsProjectile()
    {
        Vector2 target = Vector2.zero;

        GameObject projectile = GameObject.FindGameObjectWithTag("Projectile");
        if (projectile != null)
            target = new Vector2(projectile.transform.position.x, projectile.transform.position.y + 1.5f ) ;

        while (Vector2.Distance(transform.position, target) > 0.25f && projectile != null)
        {
            transform.position = Vector2.MoveTowards(transform.position, target, 1000 * Time.unscaledDeltaTime);
            yield return null;
        }
        StartCoroutine(resetIndicatorIcon());
        isWeaponInHand = true;
        canThrow = true;
        isBeganThrowing = false;
        anim.SetBool("isTeleporting", false);
        if (projectile != null)
            Destroy(GameObject.FindGameObjectWithTag("Projectile").gameObject);
        m_AttackController.canAttack = true;
        isTeleporting = false;
        if (cameraSwitch != -1)
        {
            manager.switchCameraBoundsBetweenRegion(cameraSwitch);
            cameraSwitch = -1;
        }
        if (regionSwitch != -1)
        {
            manager.setCameraBounds(regionSwitch);
            regionSwitch = -1;
        }
        
    }

    IEnumerator resetIndicatorIcon()
    {
        if(manager.tutorialBegin && manager.isFirstGame)
        {
            manager.BeginKunaiTutorial();
        }
        while (launchCoolDown > 0)
        {
            yield return null;
            manager.throwableIndicator.fillAmount = Mathf.MoveTowards(manager.throwableIndicator.fillAmount, 1, Time.deltaTime/launchRate);
            launchCoolDown -= Time.deltaTime;
        }
    }

    public void resetThrowingVars()
    {
        if(zoomout != null)
            StopCoroutine(zoomout);
        StartCoroutine(resetCamera());
        GameObject proj = GameObject.FindGameObjectWithTag("Projectile");
        if (proj != null)
        {
            StartCoroutine(resetIndicatorIcon());
            Destroy(proj);            
        }
        Time.timeScale = 1;
        Time.fixedDeltaTime = 0.02F * Time.timeScale;
       // launchCoolDown = 0;
        canThrow = true;
        retrieveWeapon = false;
        isWeaponInHand = true;
        isBeganThrowing = false;
        enemyStun = false;
        if(m_AttackController != null)
            m_AttackController.canAttack = true;
        firingHeight = minFiringHeight;
        if(dofComp != null)
            dofComp.mode.value = DepthOfFieldMode.Off;
       
        Physics2D.gravity = new Vector2(Physics2D.gravity.x, -9.8f);
        if(rb!= null)
            rb.gravityScale = 2;

        firingHeight = minFiringHeight;
        isBeganThrowing = false;

    }

    public void TeleportAndSwitchCamera(int part)
    {
        cameraSwitch = part;
    }
    public void TeleportAndSwitchRegion(int region)
    {
        regionSwitch = region;
    }

    private IEnumerator resetCamera()
    {
        
        while(cam.orthographicSize > minCamZoomSize)
        {
            yield return null;
            cam.orthographicSize = Mathf.MoveTowards(cam.orthographicSize, minCamZoomSize, 5 * Time.unscaledDeltaTime);           
        }

    }
   
}
